
#include "stdinc.h"

main(int argc, char** argv)
{
    if (argc < 3) {
	cerr << "usage: perturb ecc\n";
	exit(1);
    }

    real a = 1;
    real m = 1;
    real e = atof(argv[1]);
    real th = atof(argv[2]);

    cerr << "e = " << e << "  th = " << th << endl;

    real E = -0.5*m/a;
    real r = a*(1-e*e)/(1+e*cos(th));
    real rp = a*(1-e);
    real ra  = a*(1+e);

    cerr << "E = " << E << "  r = " << r << endl;
    cerr << "rp = " << rp << "  ra = " << ra << endl;

    real j0 = M_PI/2 + asin((r/a-1)/e);
    real y = r - rp;

    cerr << "j0 = " << j0 << "  y = " << y << endl;

    real coeff = -1.5*a*a/(e*sqrt(2*abs(E)));
    real term1 = e*e*j0;
    real term2 =  sqrt((ra-r)*(r-rp))*(2-e*e+y/(3*a));

    cerr << "coeff = " << coeff << "  term1 = " << term1
	 << "  term2 = " << term2 << endl;

    cerr << "result = " << coeff * (term1 - term2) << endl;
}
