#include <stdio.h>
#include <math.h>
#include "inverse.h"

double arc_length;

doublereal gel_(doublereal *u, doublereal *nc, doublereal *mc, doublereal *b, 
	doublereal *d__, doublereal *j, doublereal *sn, doublereal *cn, 
	doublereal *dn)
{
    /* System generated locals */
    doublereal ret_val;

    ret_val = *b + *mc * *d__ - arc_length;
    return ret_val;
} /* gel_ */

doublereal gel1_(doublereal *u, doublereal *nc, doublereal *mc, doublereal *b,
	 doublereal *d__, doublereal *j, doublereal *sn, doublereal *cn, 
	doublereal *dn, doublereal *b1, doublereal *d1, doublereal *j1, 
	doublereal *sn1, doublereal *cn1, doublereal *dn1)
{
    /* System generated locals */
    doublereal ret_val;

    ret_val = *b1 + *mc * *d1;
    return ret_val;
} /* gel1_ */

double cel(double kc, double p, double a, double b)
{
	//printf("cel() -> %f\n", cel_(&kc, &p, &a, &b));
}

double cel2(double k)
{
	double one = 1;
	double k2 = k*k;
	//printf("%f -> %f\n", k, cel_(&k, &one, &one, &k2));
	//printf("%f -> %f\n", k, cel_(&one, &one, &one, &one));
}

double perimeter(double a, double b)
{
	double one = 1;
	// find the ecentricity
	double e = sqrt(1-(b/a)*(b/a));
	// we need the complementary modulus
	double ec = sqrt(1-e*e);
	double e2 = ec*ec;
	return 4*a*cel_(&ec, &one, &one, &e2);
}

// bigel - the bisection method
// aigel - the accelerated bisection method
// nigel - Newton's method using the accelerated bisection method to prepare its starter

int main()
{
	double a = 2;
	double b = 1;
	double n = 0, m = 1-(b*b)/(a*a);
	cel(1e-1, 4.1, 1.2, 1.1);
	cel(1e-1, -4.1, 1.2, 1.1);
	cel2(1);
	double nc = 1. - n, mc = 1. - m;
	double kc=sqrt(mc);
	double one = 1., zero = 0.;
	double cB=cel_(&kc,&one,&one,&zero);
	double cD=cel_(&kc,&one,&zero,&one);
	double cJ=cel_(&kc,&nc,&zero,&one);
	double phi, u, B, D, J, f, sn, cn, dn;
	double error = 1.-15;
	//printf("%f %f -> %f\n", a, b, perimeter(a, b)/4);
	//printf("%f %f -> %f\n", a, b, perimeter(a, b)/4);
	double p = perimeter(a, b);
	double quarter_perimeter = p/4;
	arc_length = 0.5;
	u=aigel_(&nc,&mc,&cB,&cD,&cJ,&error,&error,gel_,&B,&D,&J,&sn,&cn,&dn,&f);
	phi=atan2(sn,cn);
	//printf("%12s %f %f\n", "accelerated:", phi, u);
	arc_length = 0;
        double total_length = 0;
        int i = 0;
	do {
                i++;
		total_length = i*(p/50);
		double ssn;
		double scn;
		double length = fmod(total_length, quarter_perimeter*2);
		if (length < quarter_perimeter/2) {
			arc_length = length;
			ssn = 1;
			scn = 1;
		} else if (length < quarter_perimeter) {
			arc_length = quarter_perimeter - length;
			ssn = 1;
			scn = -1;
		} else if (length < quarter_perimeter*3/2) {
			arc_length = length - quarter_perimeter;
			ssn = -1;
			scn = -1;
		} else {
			arc_length = quarter_perimeter*2 - length;
			ssn = -1;
			scn = 1;
		}
		u=nigel_(&nc,&mc,&cB,&cD,&cJ,&error,&error,gel_,gel1_,&B,&D,&J,&sn,&cn,&dn,&f);
		phi=atan2(sn,cn);
		//printf("%12s %f %f %f %f %f\n", "Newton:",arc_length, phi,u, ssn*a*sn, scn*cn);
		printf("%f %f\n", ssn*a*sn, scn*cn);
	} while (total_length*2 < p);
}


