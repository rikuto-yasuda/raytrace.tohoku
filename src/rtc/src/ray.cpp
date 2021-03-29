// tracer.cpp: tracer �N���X�̃C���v�������e�[�V����
//
//////////////////////////////////////////////////////////////////////
#include "StdAfx.h"
#include "ray.h"

using namespace rtc;

//////////////////////////////////////////////////////////////////////
// �\�z/����
//////////////////////////////////////////////////////////////////////
ray::ray( const wave_parameter& wparam )
 : m_dt_before(0.0),
   m_wave( wparam ),
   m_rk (
      boost::numeric::ublas::zero_vector<double>(3),
	  boost::numeric::ublas::zero_vector<double>(3)
   ),
   m_drk(
      boost::numeric::ublas::zero_vector<double>(3),
	  boost::numeric::ublas::zero_vector<double>(3)
   )
{}

ray::~ray()
{}

// ray::intermediate ///////////////////////////////////////////////////
ray::intermediate& ray::intermediate::operator =(
	const ray::intermediate& r
){
	B    = r.B;
	Bk   = r.Bk;
	B2   = r.B2;
	k2   = r.k2;
	w    = r.w;
	w2   = r.w2;
	wp2  = r.wp2;
	wc2  = r.wc2;
	s    = r.s;
	c    = r.c;
	s2   = r.s2;
	c2   = r.c2;
	X2   = r.X2;
	Y2   = r.Y2;
	X4   = r.X4;
	Y4   = r.Y4;
	iX2  = r.iX2;
	iX22 = r.iX22;
	numerator   = r.numerator;
	denominator = r.denominator;
	root = r.root;

	return *this;
}
bool ray::intermediate::operator ==(
	const ray::intermediate& r
) const {
	// �������̂��߁A&& �ł͂Ȃ�'&'���g���B
	// �S�� bool �� �P�ł� false ������΂����󂾂���
	// ����ł����E�E�E�͂�
	assert( B.size() == r.B.size() );
	return (
		( B(0) == r.B(0) )&
		( B(1) == r.B(1) )&
		( B(2) == r.B(2) )&
		( Bk  == r.Bk    )&
		( B2  == r.B2    )&
		( k2  == r.k2    )&
		( w   == r.w     )&
		( w2  == r.w2    )&
		( wp2 == r.wp2   )&
		( wc2 == r.wc2   )&
		( s   == r.s     )&
		( c   == r.c     )&
		( s2  == r.s2    )&
		( c2  == r.c2    )&
		( X2  == r.X2    )&
		( Y2  == r.Y2    )&
		( X4  == r.X4    )&
		( Y4  == r.Y4    )&
		( iX2 == r.iX2   )&
		( iX22 == r.iX22 )&
		( numerator   == r.numerator   )&
		( denominator == r.denominator )&
		( root == r.root )
	);
}

// �r���ŗ��p����v�Z�l�̍X�V //////////////////////////////////////////

void ray::update_intermediate(
	ray::intermediate& i,
	const vector&      r,
	const vector&      k
) const {
	const cosmos& c = getCosmos();

	i.B      = c.getMagnetField(r);
	i.Bk     = clearNaN( inner_prod(i.B,k)             );
	i.B2     = clearNaN( inner_prod(i.B,i.B)           );
	i.k2     = clearNaN( inner_prod(k,k)               );
	i.wp2    = clearNaN( c.getSquaredPlasmaFreq(r) );
	i.wc2    = clearNaN( c.getSquaredCycloFreq(r)  );

	i.w      =  m_wave.getFreq()             ;
	i.w2     =  i.w*i.w                      ;

	const double theta = std::acos(
		inner_prod( i.B/norm_2(i.B), k/norm_2(k) )
	); // theta = ����Ɛ����p
	i.s      =  std::sin(theta)              ;
	i.c      =  std::cos(theta)              ;
	i.s2     =  i.s*i.s                      ;
	i.c2     =  i.c*i.c                      ;

	i.X2     = clearNaN( i.wp2/i.w2          );
	i.Y2     = clearNaN( i.wc2/i.w2          );

	i.X4     =  i.X2*i.X2                    ;
	i.Y4     =  i.Y2*i.Y2                    ;
	i.iX2    =  1-i.X2                       ;
	i.iX22   =  i.iX2*i.iX2                  ;

	i.root        = std::sqrt(
		  i.Y4 * i.s2*i.s2
		+ 4.0  * i.Y2*i.c2*i.iX22
	);

	i.numerator   = numerator_G  (i,r)       ;
	i.denominator = denominator_G(i,r)       ;
}

// �����̏�ԃ`�F�b�N //////////////////////////////////////////////////
void ray::checkState(
	const ray::intermediate& i,
	const vector& r,
	const vector& k
) const {

	// ���������ĂȂ����ǂ����`�F�b�N����B
	if( i.numerator >= i.denominator )
		throw std::runtime_error(log("the ray was evanesced."));


	// �s����NaN�l�͗�O�G���[�𓊂���B
	const bool isNaN_r = 
		(isnan(r[0]) != 0) |
		(isnan(r[1]) != 0) |
		(isnan(r[2]) != 0) ;

	if( isNaN_r )
		throw std::runtime_error(
			log("r vector reached an invalid value.")
	);

	const bool isNaN_k =
		 (isnan(k[0]) != 0) |
		 (isnan(k[1]) != 0) |
		 (isnan(k[2]) != 0) ;

	if( isNaN_k )
		throw std::runtime_error(
			log("k vector reached an invalid value.")
	);

	// ���ܗ��� 0 < n <= 1 �ɂ��邱�Ƃ��m�F����B
	// ���Ȃ݂ɂ��͈̔͂���E�����烂�[�h�ϊ����N����炵���B
	const double n = cnst::c * norm_2(k)/( 2*cnst::pi*getWaveParam().getFreq() );
	if( n <= 0.0 || 1.0 < n )
		throw std::runtime_error(
			log("core::ray : The refractive index reached outside the range.")
	);

	//�ő̕��N�����Ɍv�Z���~�߂�
	if( (std::sqrt((pow(r(0),2.0))+(pow(r(1),2.0))+(pow(r(2)+1.601e6,2.0)))) < 1.601e6 )
		throw std::range_error(
			log("enter solid part.")
	);
}

// �g���x�N�g���̏����l ////////////////////////////////////////////////
ray* ray::initialize(
	const vector& r,
	const vector& k
) {
	// r(x,y,z)����k(x,y,z)�����������Ă���
	// �g���x�N�g���𐶐����ĕԂ��B
	intermediate i;
	update_intermediate(i,r,k);
	const vector k_new = k
		* std::sqrt( 1 - i.numerator/i.denominator )
		* m_wave.getFreq()
		/( cnst::c * norm_2(k) )
	;

	// initialize�� m_im �����������Ȃ���΂Ȃ�Ȃ��B
	update_intermediate(m_im,r,k_new);
	checkState(m_im,r,k_new);
	
	// initialize�ŁAm_dt_before�����������Ȃ���΂Ȃ�Ȃ��B
	m_dt_before = getWaveParam().getTimeStep().second;

	// k�x�N�g���̐�Βl�́A�֐�G��ό`���邱�Ƃŏo�����Ƃ��ł���B
	m_rk = vector_pair( r, k_new );
	
	return this;
}

ray* ray::initialize(
	double rx, double ry, double rz,
	double kx, double ky, double kz
) {
	assert( kx*kx + ky*ky + kz*kz != 0 );

	vector
		r = boost::numeric::ublas::zero_vector<double>(3),
		k = boost::numeric::ublas::zero_vector<double>(3);
	r(0) = rx; r(1) = ry; r(2) = rz;
	k(0) = kx; k(1) = ky; k(2) = kz;

	return initialize(r,k);
}

ray* ray::initialize(
	const vector& r,
	double pitch, double round
) {
	// r(rx,ry,rz)�ɂ����鎥��B�ɂ��āA���C���[�����g�x�N�g����N�Ƃ����
	// �x�N�g��B��(N x B)�����Ƃ���pitch�p��]���A����B������round������]����
	// ������g���x�N�g���Ƃ��Đ������A�����x�N�g���y�A��Ԃ��B
	vector
		m = getCosmos().getPlanet().getMagnet().getMagneticMoment(),
		B = getCosmos().getMagnetField(r);

	// ���ꂪ�Ȃ��Ƃ܂����̂Ń`�F�b�N�B
	double norm_m = norm_2(m), norm_B = norm_2(B);

	if( norm_m == 0.0 )
	{
		// ���ꃂ�[�����g�́A�Ƃ肠����Z�������ɐL�т���̂��g���B
		log("libraytrace: ray<>::initialize: "
			"magnetic moment not exists, instead we use vector(0,0,1).");
		m(2) = 1.0;
		norm_m = norm_2(m);
	}
	if( norm_B == 0.0 )
	{
		// ����́A�Ƃ肠����-m�������B
		B = -m;
		norm_B = norm_2(B);
	}

    // ��]�ŗ��p���邽�߂ɁA�P�ʃx�N�g��������B
    m /= norm_m;
    B /= norm_B;

	// ���i�̉�]�B���̉�]�́AN�~B�����Ƃ��āApitch[rad]������]����B
	// outer_prod()�̓e���\���ςȂ̂Ŏg���Ȃ������肷��B
	vector n = boost::numeric::ublas::zero_vector<double>(3);
	n(0) = m(1)*B(2)-m(2)*B(1),
	n(1) = m(2)*B(0)-m(0)*B(2),
	n(2) = m(0)*B(1)-m(1)*B(0);

	vector vk = rotation( B, n, pitch );
	
	// ��Q�i�̉�]�B���̉�]�́AB�����Ƃ���round[rad]������]����B
	vk = rotation( vk, B, round );

	// ��]�������ʂ��A�g���x�N�g���̕����ł���B
	return initialize(r,vk);
}

ray* ray::initialize(
	double rx, double ry, double rz,
	double pitch, double round
) {

	vector r = boost::numeric::ublas::zero_vector<double>(3);
	r(0) = rx; r(1) = ry; r(2) = rz;

	return initialize( r, pitch, round );
}

// 1step�̎��ԗʂ��v�Z /////////////////////////////////////////////////
double ray::calc_dt(
	const vector_pair&      rk,
	const vector_pair&     drk,
	ray::intermediate&   out_i
) const {
	assert( m_dt_before > 0.0 );

	// �O��x�N�g�����̋��e�͈�
	const double precision = m_wave.getPrecision();

	// step.first == max, step.second == min;
	const std::pair<double,double>& step = getWaveParam().getTimeStep();
	double        dt = m_dt_before;
	const  double abs_drdt = norm_2(drk.first);

	// dt = m_dt_before ���̑O��x�N�g�����r���A
	// ���e�͈͓��Ȃ�A���e�͈͂��肬��܂�dt���̂΂�
	// ���e�͈͊O�Ȃ�A���e�͈͂܂�dt��������B
	vector
		r = rk.first  + dt*drk.first,
		k = rk.second + dt*drk.second;

	intermediate i;
	update_intermediate(i,r,k);

	double ratio = norm_2(
		calc_dGdk(i,r,k)/calc_dGdw(i,r,k)
	) / abs_drdt;

	if(
	   1.0 - precision < ratio && ratio < 1.0 + precision &&
	   abs_drdt*dt < m_wave.getStepLength()               &&
	   dt < step.first                                    &&
	   i.numerator < i.denominator
	){
		// �덷�͈͓��ł���A����
		// �Q�{�̎��Ԃɂ��Ă��A���Ԑ������������������ł���
		// �덷�͈͓��ł����Ă��A���ԂƋ��������ɂ�����̂ł����
		// ����ȏ㎞�Ԃ𑝂₷���Ƃ͂ł��Ȃ��B
		do
		{
			out_i = i;
			dt *= 2.0;

			r = rk.first  + dt*drk.first,
			k = rk.second + dt*drk.second;
			update_intermediate(i,r,k);

			ratio = norm_2(
				calc_dGdk(i,r,k)/calc_dGdw(i,r,k)
			) / abs_drdt;

		} while (
			1.0 - precision < ratio && ratio < 1.0 + precision &&
			abs_drdt*dt < m_wave.getStepLength()               && 
			dt < step.first  /* ���ԁA�������͈͓��ł���*/     &&
			i.numerator < i.denominator
		);

#ifndef RTC_RAYTRACE_RAY_LOGS_OUT_OF_TIME_RANGE
		if( dt >= step.first )
		{ log("ray::calc_dt : out of dt range, over."); }
#endif//RTC_RAYTRACE_RAY_LOGS_OUT_OF_TIME_RANGE

		// dt �������𒴂����Ƃ���ŋA���Ă���̂�
		// dt/2�ɂ��Đ������ɖ߂��B
		dt *= 0.5;
	}

	else
	{
		assert( dt >= step.second );

		// ���������𒴂��Ă���\�������邩��
		// ���������𖞂������܂Ŏ��Ԃ����炷�B
		while(
			dt          > step.second           &&
			abs_drdt*dt > m_wave.getStepLength()
		){  dt *= 0.5; };

		// �덷�͈͊O�ł���A����
		// 1/2�{�̎��Ԃɂ��Ă��A���Ԑ������ł���B
		// �덷�͈͓��ł����Ă��A���Ԑ����ɂ�����̂ł����
		// ����ȏ㎞�Ԃ����炷���Ƃ͂ł��Ȃ��B
		// ���Ԃ����炷���\�b�h�Ȃ̂ŁA���������ɂ����邱�Ƃ͂��蓾�Ȃ��B
		out_i = i;
		do
		{
			dt *= 0.5;

			r = rk.first  + dt*drk.first,
			k = rk.second + dt*drk.second;
			update_intermediate(i,r,k);

			ratio = norm_2(
				calc_dGdk(i,r,k)/calc_dGdw(i,r,k)
			) / abs_drdt;

			if( 
			   step.second > dt               ||
			   i.numerator >= i.denominator
			){
#ifndef RTC_RAYTRACE_RAY_LOGS_OUT_OF_TIME_RANGE
				log("ray::calc_dt() : out of dt range, under.");
#endif//RTC_RAYTRACE_RAY_LOGS_OUT_OF_TIME_RANGE
				
#ifdef RTC_RAYTRACE_ENABLE_EXCEPTION_WHEN_TIMESTEP_UNDERFLOW
				throw std::runtime_error(
					"ray::calc_dt : dt range under flow."
				);
#endif//RTC_RAYTRACE_ENABLE_EXCEPTION_WHEN_TIMESTEP_UNDERFLOW
				
				// ���Ԑ�������O�ꂽ�̂ŁA2�{�ɂ��ĕԂ��B
				dt *= 2.0;
				break;
			}

			out_i = i;
		} while (
			ratio < 1.0 - precision || 1.0 + precision < ratio
		);

		// dt �����[�g�����ɓ��������A���Ԑ�������O�ꂽ�Ƃ���ŋA���Ă���
		// dt �� 0.5�{�ɂ���ƁA���[�g��������O���\��������B
	}

	assert( dt == std::min( step.first, std::max( step.second, dt ) ) );
	return m_dt_before = dt;
}


// ���߂������ ////////////////////////////////////////////////////////
double ray::calc_dGdw(
	const ray::intermediate& i,
	const vector&            r,
	const vector&            k
) const {

	const int    ox         = m_wave.LO_or_RX();
	const double first_term = -2*i.k2*(cnst::c*cnst::c)/(i.w2*i.w);

	const double
		// numerator_G()��w�Ŕ����������ʒl
		dnG_dw  = clearNaN( -4.0*(i.X2 - 2.0*i.X4)/i.w ),

		// denominator_G()��w�Ŕ����������ʒl
		ddG_dw  = clearNaN(
		          4*(i.X2/i.w)
		        + 2*(i.Y2/i.w)*i.s2
				+ ox*2*(
				  -(i.Y4/i.w) * (i.s2*i.s2)
				  -2*(i.Y2/i.w) * i.c2*(i.iX22)
				  +4*(i.X2*i.Y2/i.w) * i.c2 * i.iX2
				)/i.root
		);

	return clearNaN(
		first_term + (
		 i.denominator*dnG_dw - i.numerator*ddG_dw
		) / (i.denominator * i.denominator)
	);
}

vector ray::calc_dGdr(
	const ray::intermediate& i,
	const vector&            r,
	const vector&            k
) const {
	vector result = boost::numeric::ublas::zero_vector<double>(3);
	const cosmos& c = getCosmos();
	
	const matrix dBdr = c.getDerivativeB(r);
	const int    ox   = m_wave.LO_or_RX();
	const vector dndx = c.getDerivativeDensity(r);

	// �x�N�g�������Ȃ̂ŁA�e�������ɔ����������ʂ��i�[����B
	for( int n=0; n<3; ++n )
	{
		const double
			kdBdr =(   k(0)*dBdr(0,n) +   k(1)*dBdr(1,n) +   k(2)*dBdr(2,n) ),
			BdBdr =( i.B(0)*dBdr(0,n) + i.B(1)*dBdr(1,n) + i.B(2)*dBdr(2,n) );

		const double
			dc2 = clearNaN(
				  2.*i.Bk / (i.B2*i.B2*i.k2)
				  * (i.B2*kdBdr - i.Bk*BdBdr)
				),
			ds2  = -dc2,
			ds4  = clearNaN( -2.0 * (1.0-i.c2) * dc2 );


		const double
			dwp2w2 = clearNaN(
				(cnst::e*cnst::e) * dndx(n)/(i.w2 * cnst::e0*cnst::me)
			),
			dwc2w2 = clearNaN(
				2*(cnst::e*cnst::e)/(i.w2 * cnst::me*cnst::me) * BdBdr
			);

		const double
			first_term = clearNaN(
				(2-4*(i.X2))*(dwp2w2)
			),
			second_term = clearNaN(
				-2*(dwp2w2) - i.s2*(dwc2w2) - (i.Y2) * ds2
				  + ox* ( 
					   2*(i.s2*i.s2)*(i.Y2) * dwc2w2
					  + (i.Y4) * ds4 
					  + 4*i.c2 * i.iX22 * dwc2w2
					  + 4*(i.Y2)*(i.iX22)*dc2
					  - 8*(i.Y2)*(i.c2)*(i.iX2)*(dwp2w2)
				  ) / (2*i.root)
			);

		result(n) = clearNaN(
			( i.denominator * first_term - i.numerator * second_term )
		    / ( i.denominator * i.denominator )
		);
	}
	return result;
}

vector ray::calc_dGdk(
	const ray::intermediate& i,
	const vector&            r,
	const vector&            k
) const {
	vector result = boost::numeric::ublas::zero_vector<double>(3);

	const int    ox = m_wave.LO_or_RX();
	const double sterm_factor = clearNaN(
		i.numerator/(i.denominator*i.denominator)
	);

	for( int n=0; n<3; ++n )
	{
		const double
			dc2dk = clearNaN(
			      2*(i.Bk)/(i.B2*i.k2*i.k2)
			      *( i.B(n)*i.k2 - i.Bk*k(n) )
				),
			ds2dk = -dc2dk,
			ds4dk = -2*(1-i.c2)*dc2dk,

			first_term  = clearNaN(
				2*cnst::c*cnst::c*k(n)/i.w2
			),

			second_term = clearNaN(
				-(i.Y2)*ds2dk
				+ ox*(i.Y4*ds4dk + 4*i.Y2*i.iX22*dc2dk)/(2*i.root)
			)
		;

		result(n) = first_term - sterm_factor*second_term;
	}
	return result;
}

// =====================================================================
double ray::numerator_G(        // �֐�G�̑�Q���̕��q�̒l��Ԃ��B
	const ray::intermediate& i,
	const vector&            r
) const {
	return clearNaN( 2 * i.X2 * i.iX2 );
}

double ray::denominator_G(      // �֐�G�̑�Q���̕���̒l��Ԃ��B
	const ray::intermediate& i,
	const vector&            r
) const {
	const int ox = m_wave.LO_or_RX();

	return clearNaN(
		2*i.iX2
		- i.Y2*i.s2
		+ ox*i.root
	);
}

