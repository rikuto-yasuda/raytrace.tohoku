#include "main.h"

/////////////////////////////////////////////////////$BI91R@1$N4J0W%b%G%k$N%W%m%C%H$N$?$a0lIt=$@5(B
int plot_plasma( rtc::basic_plasma_model& p )
{
	std::clog << "plotting plasma-xy";
	plot_plasma_H("pxy_normal", p, 0 );

	std::clog << "plotting plasma-xz";
	plot_plasma_V("pxz-normal", p, 0 );
/*                                                   ////////$B%-%c%S%F%#!<$OB8:_$7$J$$$N$GL5;k!J%3%a%s%H%"%&%H!K(B
	const double Re = p.getMother().getRadius();
	p.addCavity( rtc::cavity(
		0.3,
		58,15,
		0, 6,
		7.0*Re,
        0.8*Re
	));
	std::clog << "plotting plasma-xy";
    plot_plasma_H("pxy_cav", p, 0 );

	std::clog << "plotting plasma-xz";
    plot_plasma_V("pxz-cav", p, 0 );
*/
	return 0;
}

int plot_magnet( rtc::basic_magnet_model& m )
{
	plot_magnet("m0", m, 0, 55);
	return 0;
}

int main()
{
	rtc::cosmos c( 1980,5,6, 0,0,0 );
	rtc::model::magnet::IGRF       m(4);
//	rtc::model::plasma::sato_earth p;                 ////////$B%W%i%:%^%b%G%k$r%F%9%H%b%G%k$K:9$7BX$((B
	rtc::model::plasma::test_simple p;
	rtc::planet::benchmark e(m,p);                    ////////$BOG@1%b%G%k(Bearth->benchmark$B$H:9$7BX$(!J%W%m%C%H$NC10L$,OG@1H>7B$N$?$a(B
	const double Re = e.getRadius();
	/*p.addCavity( rtc::cavity(
		0.03,
		60,10,
		0,6,
		7.0*Re,
		0.8*Re
		));*/
	c.registerPlanet(e);

	//plot_magnet(m);
	plot_plasma(p);

	return 0;
}
