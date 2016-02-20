/* inverse.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

typedef double doublereal;
typedef int integer;
//#include "f2c.h"
#include <math.h>
typedef doublereal (*D_fp)(), (*E_fp)();


doublereal gel_(doublereal *u, doublereal *nc, doublereal *mc, doublereal *b, 
	doublereal *d__, doublereal *j, doublereal *sn, doublereal *cn, 
	doublereal *dn);

doublereal gel1_(doublereal *u, doublereal *nc, doublereal *mc, doublereal *b,
	 doublereal *d__, doublereal *j, doublereal *sn, doublereal *cn, 
	doublereal *dn, doublereal *b1, doublereal *d1, doublereal *j1, 
	doublereal *sn1, doublereal *cn1, doublereal *dn1);

doublereal nigel_(doublereal *nc, doublereal *mc, doublereal *cb, doublereal *
	cd, doublereal *cj, doublereal *utol, doublereal *ftol, D_fp gel, 
	D_fp gel1, doublereal *b, doublereal *d__, doublereal *j, doublereal *
	sn, doublereal *cn, doublereal *dn, doublereal *f);

doublereal aigel_(doublereal *nc, doublereal *mc, doublereal *cb, doublereal *
	cd, doublereal *cj, doublereal *utol, doublereal *ftol, D_fp gel, 
	doublereal *b, doublereal *d__, doublereal *j, doublereal *sn, 
	doublereal *cn, doublereal *dn, doublereal *f);

double cel(double qqc, double pp, double aa, double bb);
