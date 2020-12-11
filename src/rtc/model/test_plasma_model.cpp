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
double plasma::test_simple::getDensity( const vector& point ) const       ////////////�V�����v���Y�}���f���̂��߂ɃR�����g�A�E�g
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

double plasma::test_simple::getDensity( const vector& point ) const               ///////////////�V�����v���Y�}���f���iz��������exp�Ō����j
{
	const double
		d = std::fabs(4e10*exp(-point(2)/1e4))                                  //////////////point�i�Q�j�͍��xz�im)�ł���ƍl���Ă���
		;

	return d;
}
