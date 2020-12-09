// basic_plasma_model.cpp: basic_plasma_model クラスのインプリメンテーション
//
//////////////////////////////////////////////////////////////////////
#include "StdAfx.h"
#include "test_model.h"

using namespace rtc;
using namespace model;

// test_null model ////////////////////////////////////////////////////////
double plasma::test_null_plasma::getDensity( const vector& point ) const
{
	return 0;
}

// test_simple model //////////////////////////////////////////////////////
/*
double plasma::test_simple::getDensity( const vector& point ) const       ////////////新しいプラズマモデルのためにコメントアウト
{
	const double Re = getCosmos().getPlanet().getRadius();

	vector pa = point, pb = point;
	pa(0) -= 3*Re;
	pb(0) += 3*Re;

	const double
		r = std::fabs( 
			  (point(0)*point(0) - Re*Re/2)
			+ (point(1)*point(1) - Re*Re/2)
			+ (point(2)*point(2) - Re*Re/2)
		);

	return 5e14*(
		Re/r
		+ 1.5*Re/(inner_prod(pa,pa))
		+ 1.5*Re/(inner_prod(pb,pb))
	);
}
*/

double plasma::test_simple::getDensity( const vector& point ) const               ///////////////新しいプラズマモデル（z軸方向にexpで減少）
{
	const double
		d = std::fabs(130*exp(-point(2)/10))
		;

	return d;
}
