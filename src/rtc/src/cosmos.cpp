//////////////////////////////////////////////////////////////////////
// cosmos.cpp
//
#include "StdAfx.h"
#include "cosmos.h"
using namespace rtc;

//////////////////////////////////////////////////////////////////////
// �\�z/����
//////////////////////////////////////////////////////////////////////

cosmos::cosmos(
	const int year,
	const int month,
	const int mday,
	const int hour,
	const int min,
	const int sec  // �F���̎����� UT �Ŏw�肷��B
) : m_planet       ( NULL )
{
	// ���ϐ��̃`�F�b�N
#if defined _OPENMP && defined RTC_RAYTRACE_ENABLE_LOG
	if( NULL == std::getenv("OMP_NUM_THREADS") )
	{
		// ����������OpenMP�����p����ĂȂ�
		std::clog <<
			"warning : cosmos::cosmos :\n"
			"OMP_NUM_THREADS not defined, OpenMP is not used.\n"
			"Please setenv OMP_NUM_THREADS."
		"\n";
	}
#endif
	
	// �������쐬���o�^
	{
		std::tm t;
		memset( &t, 0, sizeof(t) );
		t.tm_year = year;
		t.tm_mon  = month;
		t.tm_mday = mday;
		t.tm_hour = hour;
		t.tm_min  = min;
		t.tm_sec  = sec;

		setUniversalTime(t);
	}

	// �O���[�o���ɓo�^
	if( g_cosmo )
	{
		throw std::runtime_error(
			"cosmos::cosmos : multiple cosmos tried to be created."
		);
	}
	g_cosmo = this;
}

cosmos::~cosmos()
{
#ifdef RTC_RAYTRACE_ENABLE_LOG
	if( !m_rays.empty() )
	{
		std::clog <<
			"debug warning : cosmos::~cosmos :\n"
			"rtc::ray object not deleted by eraseRay() remains. Erase automatically."
		"\n";
	}
#endif

	// ���̔j��
	rays_t::iterator it;
	for( it = m_rays.begin(); it != m_rays.end(); ++it )
	{
		delete *it;
	}

	assert( g_cosmo );
	g_cosmo = NULL;
}


// �����̐��� //////////////////////////////////////////////////////////
ray* cosmos::createRay( const wave_parameter& wparam )
{
	ray* r = new ray(wparam);
	
	boost::mutex::scoped_lock lock( m_lock );
	m_rays.insert(r);
	
	return r;
}

ray* cosmos::createRay(
	wave_parameter::wave_mode mode,   // �g���̃��[�h�BLO_MODE��RX_MODE���w�肷��B
	double                    freq,   // �g���̎��g��[Hz]���w�肷��B
	double              prec  /* = 3.74e-4 */, // step�O��̔�̋��e�����w�肷��B
	double              lstep /* = cnst::c */, // ������1step�Ői�ލő咷���w��
	const double timeStep_max /* = 1e0     */, // 1step�Ői�ގ��Ԃ̍ő�l���w��
	const double timeStep_min /* = 1e-54   */  // 1step�Ői�ގ��Ԃ̍ŏ��l���w��
){
	return createRay( wave_parameter(
		mode,
		freq,
		prec,
		lstep,
		timeStep_max,
		timeStep_min
	));
}

void cosmos::eraseRay( ray* pray )
{
	boost::mutex::scoped_lock lock( m_lock );

	rays_t::iterator it = std::find( m_rays.begin(), m_rays.end(), pray );
	if( it != m_rays.end() )
	{
		delete *it;
		m_rays.erase( it );
	}

	else throw std::runtime_error(
		"cosmos::eraseRay : ray object not created with cosmos::createRay() tried to be erased."
	);
}

// �f���n�̓o�^ --------------------------------------------------------
bool cosmos::registerPlanet( basic_planet& planet )
{
	const bool result = ( m_planet == NULL );
	if( result ) {
		m_planet = &planet;
		planet.create();
	}
	
	return result;
}


// �����I�Ȑ��l�ւ̖₢���킹 ------------------------------------------
vector cosmos::getMagnetField( const vector& r ) const
{ return m_planet->getMagnet()(r); }

// getDerivativeB() -----------------------------
// �h�����ꂽ�ł��낤 getField()���Ăяo���A
// 0.5[m]�����O��Ɉړ����Ď���̍������Ƃ߁A
// ������z�𓱏o����B
matrix cosmos::getDerivativeB( const vector& r ) const
{
	matrix dbdr(3,3);

	for( int n=0; n<3; ++n )
	{
		vector dr = boost::numeric::ublas::zero_vector<double>(3);
		dr(n)=RTC_DERIVATIVE_DISTANCE;

		const vector dB = getMagnetField(r+dr)-getMagnetField(r-dr);

		for( int i=0; i<3; ++i )// i = x,y,z
		{
			dbdr(i,n) = dB(i);
		}
	}
	return dbdr;
}

// =====================================================================
vector cosmos::getDerivativeDensity( const vector& r ) const
{
	vector result = boost::numeric::ublas::zero_vector<double>(3);

	for( int n=0; n<3; ++n )
	{
		vector ra = r, rb = r;
		ra(n) += RTC_DERIVATIVE_DISTANCE/2.0;
		rb(n) -= RTC_DERIVATIVE_DISTANCE/2.0;

		result(n) = clearNaN(
			(getCosmos().getPlasmaDensity(ra) - getCosmos().getPlasmaDensity(rb))
			/ RTC_DERIVATIVE_DISTANCE
		);
	}
	return result;
}


double cosmos::getPlasmaDensity( const vector& r ) const
{ return m_planet->getPlasma()(r); }

double cosmos::getHight( const vector& r ) const   //////���˗p�@�\�@�W������́i�ő́E�v���Y�}���E�ʂ���̋����𓱏o�j�@�P�ʂ�m  
{ 
	double h = r(2)-1.0e4;  ///z=0�̕���
///	double h = std::sqrt((pow(r(0),2.0))+(pow(r(1),2.0))+(pow(r(2)+1.601e6,2.0))) - 1.601e6;�@// �G�E���p�\�ʁi���S���j
	return h; 
}


// ���f������v�Z���Č��ʂ�Ԃ��₢���킹 //////////////////////////////
double cosmos::getSquaredPlasmaFreq( const vector& r ) const
{
	return getPlasmaDensity(r)
		* (cnst::e*cnst::e)
		/ (cnst::e0*cnst::me);
}

double cosmos::getSquaredCycloFreq ( const vector& r ) const
{
	const vector B  = getMagnetField(r);
	const double B2 = inner_prod(B,B);

	return B2
		* (cnst::e*cnst::e)
		/ (cnst::me*cnst::me);
}

