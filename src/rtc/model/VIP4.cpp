////////////////////////////////////////////////////////////////////////
// VIP4.cpp
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
#include "StdAfx.h"
#include "VIP4.h"
#include "vip4_coeffs.h"

using namespace rtc;
using namespace rtc::model::magnet;

// VIP4 - master =====================================================
VIP4::VIP4( const int dimension )
: m_coefficients( dimension )
{}

int VIP4::create( basic_planet& mother )
{
	// ��V�̂̎��Ɉʒu���w�肷��B
	const double
		lng = std::atan2( m_coefficients.h(1,1), m_coefficients.g(1,1) ),
		lat = 0.5*cnst::pi - std::atan2(
			m_coefficients.g(1,1) * std::cos(lng)
			+ m_coefficients.h(1,1) * std::sin(lng),
			m_coefficients.g(1,0)
		);
	mother.setMagneticalAxis(
		rad2deg(lat), rad2deg(lng)
	);

	// ��]�s����쐬����B
	m_sm2geo = boost::numeric::ublas::prod(
		mother.getGEI2GEO(),
		makeMatrixInverse( mother.getGEI2GSE() )
	);
	m_sm2geo = boost::numeric::ublas::prod(
		m_sm2geo,
		makeMatrixInverse( mother.getGSE2GSM() )
	);
	m_sm2geo = boost::numeric::ublas::prod(
		m_sm2geo,
		makeMatrixInverse( mother.getGSM2SM() )
	);
	m_geo2sm = makeMatrixInverse( m_sm2geo );

	return basic_magnet_model::create( mother );
}

vector VIP4::getField( const vector& pos ) const
{
	const basic_planet& earth = getMother();
	vector r = boost::numeric::ublas::zero_vector<double>(3);

	// pos�ɂ� SM�������W�n�ł̈ʒu�x�N�g�����n�����̂Œ��ӂ���B
	// VIP4�ł́AGEO�����W�n�ɕϊ����Ȃ���΂Ȃ�Ȃ��B
	// �܂��A�߂�l��SM�������W�n�Ŗ�����΂Ȃ�Ȃ��B
	const vector ptr = convertToPolar(
		rotation_prod( m_sm2geo, pos )
	);
#ifdef RTC_RAYTRACE_ENABLE_DETAIL_LOG
	std::clog << "VIP4 convert: from " << pos << "\n"
			  << "to " << 90.0-rad2deg(ptr[1]) << ", " << rad2deg(ptr[2]) << "\n";
#endif
	
	const double 
		c_t = std::cos(ptr[1]),
		s_t = std::sin(ptr[1]),
		ratio = earth.getRadius() / ptr[0];

	const int
		n_max = m_coefficients.getDimensionEnd();


	boost::multi_array<double,2>
		P( boost::extents[n_max][n_max] ),
		Q( boost::extents[n_max][n_max] );
	P[0][0] = 1.0;
	P[1][0] = c_t;
	Q[0][0] = 0.0;

	// ���ɂ̐^��ɂ���ꍇ�As_t == 0.0�ƂȂ�O���Z�G���[����������B
	if( s_t != 0.0 )
	{
		for( int n = 1 ;; )
		{
			const double
				powerd_r = std::pow( ratio, n+2.0 );

			Q[n][0] = n*(c_t*P[n][0] - P[n-1][0])/s_t;
			
			for( int m = 0 ;; )
			{
				const double
					c_ml      = std::cos( m * ptr[2] ),
					s_ml      = std::sin( m * ptr[2] ),
					G        = m_coefficients.G(n,m),
					H        = m_coefficients.H(n,m);

				r[0] += powerd_r * ( G*c_ml + H*s_ml )         * Q[n][m]; //  X
				r[1] += powerd_r * m * ( G*s_ml - H*c_ml )     * P[n][m]; //  Y
				r[2] -= powerd_r * (n+1) * ( G*c_ml + H*s_ml ) * P[n][m]; // -Z
				
				// �I������
				if( ++m <= n )
				{
					P[n][m] = (
					    (n+m-1)*P[n-1][m-1] - (n-m+1)*c_t*P[n][m-1]
				    )/s_t;
					Q[n][m] = (n+m)*(n-m+1) * P[n][m-1] - m*c_t/s_t * P[n][m];
					continue;
				}
				else break;
			}
			
			if( ++n < n_max )
			{
				P[n][0] = (
			   		(2*n - 1)*c_t*P[n-1][0] - (n-1)*P[n-2][0]
				) / n;
				
				continue;
			}
			else break;
		}
		r[1] /= s_t;
	}

	// ���ɂ̐^��ɂ���ꍇ�AP == 1.0, Q == 0.0�ł���B
	else
	{
		for( int n = 1; n < n_max; ++n )
		{
			const double
				powerd_r = std::pow( ratio, n+2.0 );

			for( int m = 0 ; m <= n; ++m )
			{
				const double
					c_ml      = std::cos( m * ptr[2] ),
					s_ml      = std::sin( m * ptr[2] ),
					G        = m_coefficients.G(n,m),
					H        = m_coefficients.H(n,m);

				// r[0] += 0.0; //  X
				r[1] += powerd_r * m * ( G*s_ml - H*c_ml )    ; //  Y
				r[2] -= powerd_r * (n+1) * ( G*c_ml + H*s_ml ); // -Z
			}
		}
	}


	// �P�ʂ� [nT] �i�����j�Ȃ̂ŁA[T]�ɕϊ�����B
	r *= 1e-9;

	// ���ʂ́AX=�k�AY=�o�x�AZ=�n�������̎���Ȃ̂�
	// ��]���� GEO���W�n�ɂ���B
	// ���̉�]�́A�܂�Y�����S�� ��-�� ������]���A
	// ����Z�����S�� �� ������]����B
	//
	// ���̌��SM���W�n�ɕϊ����A�Ԃ��B

	matrix rm = boost::numeric::ublas::prod(
		m_geo2sm,
		makeMatrixRotationZ( ptr[2] )
	);
	rm = boost::numeric::ublas::prod(
		rm,
		makeMatrixRotationY( ptr[1] - cnst::pi )
	);

	return rotation_prod( rm, r );
}


// VIP4 - coefficient ================================================
VIP4::coefficient::coefficient( const int n_max )
: m_elements( boost::extents[ coeff_array_t::extent_range(1,n_max+1)][n_max+1] )
{
	load();
}

double VIP4::coefficient::g( int n, int m ) const
{
	assert( 1 <= n && n < getDimensionEnd() );
	assert( 0 <= m && m <= n );
	return m_elements[n][m].g;
}

double VIP4::coefficient::h( int n, int m ) const
{
	assert( 1 <= n && n < getDimensionEnd() );
	assert( 0 <= m && m <= n );
	return m_elements[n][m].h;
}

double VIP4::coefficient::G( int n, int m ) const
{
	assert( 1 <= n && n < getDimensionEnd() );
	assert( 0 <= m && m <= n );
	return m_elements[n][m].G;
}

double VIP4::coefficient::H( int n, int m ) const
{
	assert( 1 <= n && n < getDimensionEnd() );
	assert( 0 <= m && m <= n );
	return m_elements[n][m].H;
}

void VIP4::coefficient::load()
{
	int c = 0;// �J�E���^��p�ϐ�

	const int max_n = getDimensionEnd();

	for( c = 0; c < g_vip4_coeffs_nums; ++c )
	{
		const int n = g_vip4_coeffs[c].n;
		const int m = g_vip4_coeffs[c].m;

		if( n >= getDimensionEnd() ){
			continue;
		}

		//  G = g[n][m] * sqrt( em * (n-m)! / (n+m)! );
		//  H = h[n][m] * sqrt( em * (n-m)! / (n+m)! );
		//  em = ( m == 0 ? 1 : 2 )
		// �̒l���i�[����B
		int i;
		assert( n >= m );
		const double em = ( m == 0 ? 1 : 2 );
		unsigned long
			nm1 = 1, // (n-m)!
			nm2 = 1; // (n+m)!
		for( i = n+m; i > n-m; --i ) {
			nm2 *= i;
		}
		for( ; i > 0; --i ) {
			nm1 *= i;
			nm2 *= i;
		}
		const double factor = std::sqrt(
			em * nm1/nm2
		);

		switch( g_vip4_coeffs[c].mode )
		{
			case VIP4_COEFFICIENT_ELEMENT::g:
				assert( m_elements[n][m].g == 0.0 );
				assert( m_elements[n][m].G == 0.0 );
				m_elements[n][m].g = g_vip4_coeffs[c].data;
				m_elements[n][m].G = m_elements[n][m].g * factor;
				break;

			case VIP4_COEFFICIENT_ELEMENT::h:
				assert( m_elements[n][m].h == 0.0 );
				assert( m_elements[n][m].H == 0.0 );
				m_elements[n][m].h = g_vip4_coeffs[c].data;
				m_elements[n][m].H = m_elements[n][m].h * factor;
				break;
		}
	}
}
