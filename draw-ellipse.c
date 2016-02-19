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
	double a = 3;
	double b = 2;
	double n = 0, m = 1-(b*b)/(a*a);
	double nc = 1. - n, mc = 1. - m;
	double kc=sqrt(mc);
	double one = 1., zero = 0.;
	double cB=cel_(&kc,&one,&one,&zero);
	double cD=cel_(&kc,&one,&zero,&one);
	double cJ=cel_(&kc,&nc,&zero,&one);
	double u, B, D, J, f, sn, cn, dn;
	double error = 1.-15;
	double p = perimeter(a, b);
	double quarter_perimeter = p/4;
        double total_length = 0;
        int i = 0;
        int number_of_points = 35;
	do {
                i++;
		total_length = i*(p/number_of_points);
		double ssn;
		double scn;
		double length = fmod(total_length, p);
		if (length < quarter_perimeter) {
			arc_length = length/a;
			ssn = 1;
			scn = 1;
		} else if (length < quarter_perimeter*2) {
			arc_length = (quarter_perimeter*2 - length)/a;
			ssn = 1;
			scn = -1;
		} else if (length < quarter_perimeter*3) {
			arc_length = (length - quarter_perimeter*2)/a;
			ssn = -1;
			scn = -1;
		} else {
			arc_length = (p - length)/a;
			ssn = -1;
			scn = 1;
		}
		u=nigel_(&nc,&mc,&cB,&cD,&cJ,&error,&error,gel_,gel1_,&B,&D,&J,&sn,&cn,&dn,&f);
		//double phi=atan2(sn,cn);
		//printf("%12s %f %f %f %f %f\n", "Newton:",arc_length, phi,u, ssn*a*sn, scn*cn);
		printf("%f %f\n", ssn*a*sn, scn*b*cn);
	} while (total_length < p);
}


