#include <stdio.h>
#include <math.h>
#include "inverse.h"
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
	u=aigel_(&nc,&mc,&cB,&cD,&cJ,&error,&error,gel_,&B,&D,&J,&sn,&cn,&dn,&f);
	phi=atan2(sn,cn);
	printf("%12s %f %f\n", "accelerated:", phi, u);
	u=nigel_(&nc,&mc,&cB,&cD,&cJ,&error,&error,gel_,gel1_,&B,&D,&J,&sn,&cn,&dn,&f);
	phi=atan2(sn,cn);
	printf("%12s %f %f\n", "Newton:",phi,u);
}


