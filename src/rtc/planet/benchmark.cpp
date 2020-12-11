////////////////////////////////////////////////////////////////////////
// benchmark.cpp
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
#include "StdAfx.h"
#include "benchmark.h"

using namespace rtc;
using namespace rtc::planet;

benchmark::benchmark(
	basic_magnet_model&   mag,        // ���ꃂ�f���̃C���X�^���X���w�肷��B
	basic_plasma_model&  plsm         // �v���Y�}���f���̃C���X�^���X���w�肷��B
) : basic_planet (
	1e3, // �n�����a[m]
//	8.43e10,  // �n���C�̉��z�o�Ɏq���[�����g[Am^2]�i�k�����Ő��j
	1e20,  // �n���C�̉��z�o�Ɏq���[�����g[Am^2]�i�k�����Ő��j
	axis_info( 90.0, 0.0 ), // �����̈ʒu
	mag,
	plsm
)
{};

/*
matrix benchmark::getGEI2GEO() const
{
	// Z���𒆐S�ɉ�]����s���n���B
	// ��]�p�����Ƃ߂鎮�́Areference�̒����Q�ƁB
	// reference�̒��ł́AT1�Ƃ��Ē�`����Ă���s��ł���B

	const std::tm& t = getCosmos().getUniversalTime();
	const double
		T0 = (getMJD() - 51544.5)/ 36525.0,
		H  = t.tm_hour + (t.tm_min/60.0) + (t.tm_sec/3600);

	const double theta = deg2rad(
		100.46
		+ 36000.770 * T0
		+ 15.04107  * H
	);

	return makeMatrixRotationZ(-theta);
}

matrix benchmark::getGEI2GSE() const
{
	// X���𒆐S�� epsilon ��]������ɁA
	// Z���𒆐S�� lambda ��]����B
	// reference�̒��ł́AT2�Ƃ��Ē�`����Ă���s��ł���B
	const std::tm& t = getCosmos().getUniversalTime();

	const double
		T0 = (getMJD() - 51544.5)/ 36525.0,
		H  = t.tm_hour + (t.tm_min/60.0) + (t.tm_sec/3600);

	const double
		M = 357.528 + 35999.050*T0 + 0.04107*H,
		A = 280.460 + 36000.772*T0 + 0.04107*H;

	const double
		epsilon = 23.439 - 0.013*T0,
		lambda  = A + (1.915 - 0.0048*T0)*std::sin( deg2rad(M) ) + 0.020*std::sin( deg2rad(2*M) );

	return boost::numeric::ublas::prod( 
		makeMatrixRotationZ( deg2rad(-lambda)  ),
		makeMatrixRotationX( deg2rad(-epsilon) )
	);
}
*/
matrix benchmark::getGSE2GSM() const
{
	// X���𒆐S�ɁApsi������]����s���Ԃ��B
	// psi �� GSE���W�ł̎����������狁�߂�B
	const vector re = getMagneticalAxisInGSE();
	assert( re[2] >= 0.0 );

	const double psi = std::atan2( re[1], re[2] );
	return makeMatrixRotationX(psi);
}

matrix benchmark::getGSM2SM () const
{
	// Y���𒆐S�ɁA-mu������]����s���Ԃ��B
	// mu �� GSE���W�n�ł̎����������狁�߂�B
	const vector re = getMagneticalAxisInGSE();
	assert( re[2] >= 0.0 );

	const double mu = std::atan2(
		re[0],
		std::sqrt( re[1]*re[1] + re[2]*re[2] )
	);

	return makeMatrixRotationY(-mu);
}

