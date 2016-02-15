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


int main()
{
	double        nc=1., mc=0.5;
	double kc=sqrt(mc);
	double one = 1., zero = 0.;
	double cB=cel_(&kc,&one,&one,&zero);
	double cD=cel_(&kc,&one,&zero,&one);
	double cJ=cel_(&kc,&nc,&zero,&one);
	double phi, u, B, D, J, f, sn, cn, dn;
	double error = 1.-15;
	arc_length = 1;
	u=aigel_(&nc,&mc,&cB,&cD,&cJ,&error,&error,gel_,&B,&D,&J,&sn,&cn,&dn,&f);
	phi=atan2(sn,cn);
	printf("%12s %f %f\n", "accelerated:", phi, u);
	arc_length = 0;
	for (int i = 0; i < 150; i++) {
		arc_length = i*9.2/149;
		u=nigel_(&nc,&mc,&cB,&cD,&cJ,&error,&error,gel_,gel1_,&B,&D,&J,&sn,&cn,&dn,&f);
		phi=atan2(sn,cn);
		printf("%12s %f %f %f %f %f\n", "Newton:",arc_length, phi,u, 2*sn, cn);
	}
}


