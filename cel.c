/* cel.f -- translated by f2c (version 20100827).
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


/*     Returns the general complete elliptic integral */

doublereal cel_(doublereal *qqc, doublereal *pp, doublereal *aa, doublereal *
	bb)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    ///* Subroutine */ int s_paus(char *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static doublereal a, b, e, f, g, p, q, em, qc;

    if (*qqc == 0.f) {
//	s_paus("failure in CEL", (ftnlen)14);
    }
    qc = fabs(*qqc);
    a = *aa;
    b = *bb;
    p = *pp;
    e = qc;
    em = 1.f;
    if (p > 0.f) {
	p = sqrt(p);
	b /= p;
    } else {
	f = qc * qc;
	q = 1.f - f;
	g = 1.f - p;
	f -= p;
	q *= b - a * p;
	p = sqrt(f / g);
	a = (a - b) / g;
	b = -q / (g * g * p) + a * p;
    }
L1:
    f = a;
    a += b / p;
    g = e / p;
    b += f * g;
    b += b;
    p = g + p;
    g = em;
    em = qc + em;
    if ((d__1 = g - qc, fabs(d__1)) > g * 3e-4f) {
	qc = sqrt(e);
	qc += qc;
	e = qc * em;
	goto L1;
    }
    ret_val = (b + a * em) * 1.5707963268f / (em * (em + p));
    return ret_val;
} /* cel_ */

