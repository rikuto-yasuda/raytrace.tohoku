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
		r = std::sqrt((pow(point(0),2.0))+(pow(point(1),2.0))+(pow(point(2)+3.2e6,2.0)));
	const double
		rxy = std::sqrt((pow(point(0),2.0))+(pow(point(1),2.0)));
	const double
		plume = std::fabs(1.0e12*exp(-(r-3.2e6)/1.5e5)*exp(-((atan2(rxy,point(2)))/0.261799)*((atan2(rxy,point(2)))/0.261799)));
//////////////エウロパ静水圧平行モデル 地表面で9.0*10^3(/cc) スケールハイト240km + プルーム
//>->---d = std::fabs(4e10*exp(-point(2)/1e4))                                  //////////////point（２）は高度z（m)・単位は（/m＾３）であると考えている
	const double
		t = std::fabs(9e9*exp(-point(2)/2.4e5));                                  //////////////エウロパ静水圧平行モデル 地表面で9.0*10^3(/cc) スケールハイト240km
//		t = std::fabs(0.5e5*point(2))                                        //////////////解析解比較用プラズマモデル
	const double
		d = t+plume;
		;

	return d;
}
