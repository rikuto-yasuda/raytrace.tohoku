////////////////////////////////////////////////////////////////////////
// basic_planet.cpp
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
#include "StdAfx.h"
#include "plasma_cavity.h"

using namespace rtc;

// cavity //////////////////////////////////////////////////////////////
cavity::cavity()
 : m_ilat_center_ptr(-1),
   m_mlt_center_ptr (-1),
   m_ilat_range     (0),
   m_mlt_range      (0),
   m_min_fpfc       (-1),
   m_max_height     (0),
   m_bottomHeight   (0),
   m_scale          (1.0),
   m_mother         (NULL)
{}

cavity::cavity(
	double min_fpfc,
	double ilat_center,
	double ilat_range,
	double mlt_center,
	double mlt_range,
	double max_height,
	double bottom_height
) : m_ilat_center_ptr(ilat_center),
    m_mlt_center_ptr (mlt_center ),
    m_ilat_range     (ilat_range ),
    m_mlt_range      (mlt_range  ),
    m_min_fpfc       (min_fpfc   ),
    m_max_height     (max_height ),
    m_bottomHeight   (bottom_height),
    m_scale          (       -1.0),
    m_mother         (       NULL)
{}

cavity& cavity::operator =(
	const cavity& c
){
	m_min_fpfc        = c.m_min_fpfc        ;
	m_ilat_center_ptr = c.m_ilat_center_ptr ;
	m_mlt_center_ptr  = c.m_mlt_center_ptr  ;
	m_ilat_range      = c.m_ilat_range      ;
	m_mlt_range       = c.m_mlt_range       ;
	m_max_height      = c.m_max_height      ;
	m_bottomHeight    = c.m_bottomHeight    ;
	m_scale           = c.m_scale           ;
	m_mother          = c.m_mother          ;

	return *this;
}

double cavity::factor( const vector& ptr ) const
{
	assert( m_mother );
	
	assert( 90 > m_ilat_center_ptr && m_ilat_center_ptr > 0 );
	assert( 24 > m_mlt_center_ptr  && m_mlt_center_ptr >= 0 );
	assert( m_ilat_range      > 0 );
	assert( m_mlt_range       > 0 );

	const double
		Re       = m_mother->getRadius(),
		ilat_max = m_ilat_center_ptr + m_ilat_range,
		ilat_min = m_ilat_center_ptr - m_ilat_range,
		mlt_max  =  m_mlt_center_ptr +  m_mlt_range,
		mlt_min  =  m_mlt_center_ptr -  m_mlt_range;
	double
		h        = m_mother->restToSoil(ptr),
		ilat     = m_mother->getFLAT(ptr, 1e-3*Re),
		mlt      = m_mother->getMLT(ptr);
	
	// �ȉ���goto���́Acavity�����ł̌Ăяo�����������낤
	// �Ƃ����\���̌��̍�������_���āB

	// m_scale�����̎��́A�L���r�e�B�������̈�ł���B
	// create()���܂��Ăяo����Ă��Ȃ����A
	// �����ȃL���r�e�B�E�p�����[�^���w�肳��Ă���Ƃ�
	// �i�R���X�g���N�^�̒l���̂܂܂Ƃ��j�̏ꍇ�A���̒l�ɂȂ�B
	if( m_scale < 0 ) {
		goto outof_cavity;
	}

	// �܂��AILAT���͈͊O�̎��� 1.0��Ԃ��B
	if( ilat < ilat_min || ilat_max < ilat ) {
		goto outof_cavity;
	}
	
	// ���ɍ���������B
	if( m_max_height < h ) {
		goto outof_cavity;
	}
	
	// �Ō��MLT������B
	if( mlt_min <  0.0 && mlt_max < mlt ){ mlt -= 24; }
	if( mlt_max > 24.0 && mlt_min > mlt ){ mlt += 24; }
	if( mlt < mlt_min  || mlt_max < mlt ){
		goto outof_cavity;
	}

	// �����܂ŗ�����A�w��ʒu��cavity�����ɂ��邱�ƂɂȂ�B
	// cavity�̌W���Ƃ��āA�����ł͎��̊֐����g�p����B
	// m_bottomHeight = A�Am_maxHeight=B�Am_scale=S�ƒu���ƁA
	//
	// case: 0.0 < h < A
	//    s(h) = cos( h*PI/(2*B) )
	//
	// case: A < h < B
	//    s(h) = sin( (h-A)*PI/(2*(B-A)) )
	//
	// ������g���āA
	//    f(h) = (1-S)*s(h)^2 + S
	// �Ƃ���B����͍��������Ɋ��炩�Ȋ֐��ł���B
	//
	// ���ɁA
	//    g(h,x,y) = (1-f(h))*( sin(x)^2 * cos(y)^2 + sin(y)^2 ) + f(h);
	// ���v�Z���A���̒l���W���Ƃ��ĕԂ��B
	// sin(x)^2 * cos(y)^2 + sin(y)^2 �́A���炩�Ȃ��o�^�̊֐��ł���B
	//
	// ���ʓ���ꂽg(h,x,y)�́A(h,x,y)��cavity�̒��S�̎��Œ�lS(0 < S < 1.0)���Ƃ�
	// ���炩�Ȋ֐��ł���B���Ȃ݂ɕ����w�I�ȈӖ��͂���܂�Ȃ��B
	// �ł� calvert et al.[1981] �Ɏ������ʂ�������B
	//
	{
		assert( h < m_max_height );
		assert( 0.0 <= m_scale && m_scale <= 1.0 );
		
		const double h_factor = (
			h < m_bottomHeight
		) ?
			std::cos( h * cnst::pi / (2*m_bottomHeight) ) // 0 < h < bottom�Ȏ�
		  :
			std::sin(
			    ( ( h - m_bottomHeight)*cnst::pi  )
				/ ( 2*(m_max_height-m_bottomHeight) )
			) // bottom <= h <= max�Ȏ�
		;
		
		//���������ɂȂ߂炩�Ȋ֐�
		const double base_scale = (1-m_scale)*(h_factor*h_factor)+m_scale;

		const double
			x = cnst::pi * (ilat-m_ilat_center_ptr)/(2*m_ilat_range),
			y = cnst::pi * ( mlt- m_mlt_center_ptr)/( 2*m_mlt_range),
			sx = std::sin(x),
			cy = std::cos(y),
			sy = std::sin(y),
			factor = (1.0-base_scale)*(
				(sx*sx * cy*cy + sy*sy) //MLT,ILAT�����Ɋ��炩�Ȋ֐�
			) + base_scale;

		// ���l��factor�Ƃ��ĕԂ��B
		return factor*factor;
	}
	
outof_cavity:
	return 1.0;
}

void cavity::create(
	 const basic_planet& mother
){
	// ���ꃂ�f�����g���[�X���āAfp/fc < 0.03 �ɂȂ�悤�ȌW��������
	// m_scale > 0 �̎��́A����create()���Ăяo����Ă���B
	if( m_scale > 0 ) {
		return;
	}
	assert( m_mother == NULL );
	m_mother = &mother;
	
	
	const double Re = mother.getRadius();
	const vector bottom_ptr = mother.getLocation(
		m_ilat_center_ptr,
		m_mlt_center_ptr,
		m_bottomHeight,
		1e-2*Re,
		NULL
	);

	// bottom_height�̓_�� fp/fc���Œ�ƂȂ�X�P�[���l���v�Z
	// fp/fc �� �w��l�ȉ��ɂȂ�悤�ɃX�P�[����ύX������
	const double fpfc = 
		getCosmos().getSquaredPlasmaFreq (bottom_ptr) 
		/ getCosmos().getSquaredCycloFreq(bottom_ptr)
	;

	// �Œ�l���Ƃ�n�_��fp/fc�ƁA����ł̒l�̔���Ƃ�
	// ���̓��l��scale��ݒ肷��B
	const double scl = m_min_fpfc/fpfc;
	m_scale = scl*scl;

#ifdef RTC_RAYTRACE_ENABLE_LOG
	std::clog <<
		"plasma cavity scale is : " << m_scale    << "\n"
		"      current fp/fc is : " << fpfc       << "\n"
		"        and minimum is : " << std::sqrt(
			getCosmos().getSquaredPlasmaFreq (bottom_ptr) 
			/ getCosmos().getSquaredCycloFreq(bottom_ptr)
		) << std::endl;
#endif
}
