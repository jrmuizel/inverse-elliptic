#include <math.h>


/*     Returns the general complete elliptic integral */

double cel(double qqc, double pp, double aa, double bb)
{
    /* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    double a, b, e, f, g, p, q, em, qc;

    if (qqc == 0.f) {
//	s_paus("failure in CEL", (ftnlen)14);
    }
    qc = fabs(qqc);
    a = aa;
    b = bb;
    p = pp;
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

