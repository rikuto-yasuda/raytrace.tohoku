////////////////////////////////////////////////////////////////////////
// nsumei_earth.cpp
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
#include "StdAfx.h"
#include "model_nsumei_earth.h"

using namespace rtc;
using namespace model;

plasma::nsumei_earth::nsumei_earth(
	const double kp
) : m_Kp   (kp)
{}

double plasma::nsumei_earth::getDensity( const vector& point ) const
{
	const double ILAT = getILAT(point);

	// �P�ʂ����킹�Ȃ���΂Ȃ�Ȃ��B
	// ���̎��ł́AR �� �n�����a�P�ʁANe ��[�R/cc]�ƂȂ��Ă���B
	const double R = norm_2(point) / getCosmos().getPlanet().getRadius();

	if( 1.4 < R && R < 5.0 &&
	    70 < ILAT
	){
		const double
			N0    =  68510,
			alpha = -0.038,
			beta  =  0.220,
			gamma =  4.95;
		return 1e6 *(
			 N0 * std::pow(R, -gamma) * std::exp( beta*m_Kp + alpha*ILAT)
		);
	}

	else {
		const double
			N0    = 3433.0,
			beta  = 0.229,
			gamma = 5.09;

		return 1e6 * (
			N0 * std::pow(R, -gamma) * std::exp( beta * m_Kp )
		);
	}
}
double plasma::nsumei_earth::getILAT( const vector& point ) const
{
	// ���[�N���b�h��� point ����A�s�ώ��C�ܓx(ILAT)�����߁A�Ԃ��B
	const vector fp = convertToPolar(
		getCosmos().getPlanet().getFootPrint(
			point,
			3e-2*getCosmos().getPlanet().getRadius()
		)
	);

	const double ilat = std::fabs( rad2deg( fp(1) )-90.0 );

	// --fix me--
	// ilat �� 90.0000000001���̐��l�ɂȂ��Ă��邱�Ƃ�����B
	// ����ɂ��e��������邽�߁A90.0�ȏ�̏ꍇ��
	// ���ׂ� 89.9999999999999 �ɂ��ĕԂ��B
	// FPU�ɂ���ẮA����ɓ��삵�Ȃ��Ȃ鋰�������̂Œ��ӂ��邱�ƁB
	return std::min( ilat, 89.9999999999999 );
}

