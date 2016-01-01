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
/* Common Block Declarations */

struct {
    integer icase;
} icase_;

#define icase_1 icase_

/* Table of constant values */

static doublereal c_b4 = .2;
static integer c__9 = 9;
static integer c__1 = 1;
static doublereal c_b15 = 0.;
static doublereal c_b20 = 1.;

doublereal gel_(doublereal *u, doublereal *nc, doublereal *mc, doublereal *b, 
	doublereal *d__, doublereal *j, doublereal *sn, doublereal *cn, 
	doublereal *dn)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal m, n;

    n = 1. - *nc;
    m = 1. - *mc;
    ret_val = *b + *mc * *d__ - 1.;
    return ret_val;
} /* gel_ */

doublereal gel1_(doublereal *u, doublereal *nc, doublereal *mc, doublereal *b,
	 doublereal *d__, doublereal *j, doublereal *sn, doublereal *cn, 
	doublereal *dn, doublereal *b1, doublereal *d1, doublereal *j1, 
	doublereal *sn1, doublereal *cn1, doublereal *dn1)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal m, n;

    n = 1. - *nc;
    m = 1. - *mc;
    ret_val = *b1 + *mc * *d1;
    return ret_val;
} /* gel1_ */

doublereal nigel_(doublereal *nc, doublereal *mc, doublereal *cb, doublereal *
	cd, doublereal *cj, doublereal *utol, doublereal *ftol, D_fp gel, 
	D_fp gel1, doublereal *b, doublereal *d__, doublereal *j, doublereal *
	sn, doublereal *cn, doublereal *dn, doublereal *f)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal);
    //integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
    //	    e_wsle(void);

    /* Local variables */
    static doublereal h__;
    static integer i__;
    static doublereal m, n, t, u, v, w, x, y, z__, b1, d1, f1, j1, bv, dv, jv,
	     xi, nu, xv, yv, zv, cn1, dn1, sn1, eta, cnv, dnv, snv, zeta;
    extern doublereal aigel_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, D_fp, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), uatan_(doublereal *, 
	    doublereal *);
    extern /* Subroutine */ int sersdj_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    //static cilist io___36 = { 0, 6, 0, 0, 0 };


    m = 1. - *mc;
    n = 1. - *nc;
    h__ = n * *nc * (*mc - *nc);
    u = aigel_(nc, mc, cb, cd, cj, &c_b4, ftol, (D_fp)gel, b, d__, j, sn, cn, 
	    dn, f);
    x = *cn * *cn;
    y = *sn * *sn;
    z__ = *dn * *dn;
    for (i__ = 1; i__ <= 10; ++i__) {
	d1 = y;
	b1 = x;
	j1 = y / (*nc + n * x);
	cn1 = -(*sn) * *dn;
	sn1 = *cn * *dn;
	dn1 = -m * *sn * *cn;
	f1 = (*gel1)(&u, nc, mc, b, d__, j, sn, cn, dn, &b1, &d1, &j1, &sn1, &
		cn1, &dn1);
	v = -(*f) / f1;
	u += v;
	sersdj_(&v, &n, &m, &snv, &dv, &jv);
	yv = snv * snv;
	xv = 1. - yv;
	zv = 1. - m * yv;
	cnv = sqrt(xv);
	dnv = sqrt(zv);
	bv = v - dv;
	xi = *cn * cnv;
	eta = *sn * snv;
	zeta = *dn * dnv;
	nu = 1. / (1. - m * y * yv);
	*cn = (xi - eta * zeta) * nu;
	*sn = (*sn * cnv * dnv + snv * sn1) * nu;
	*dn = (zeta - m * eta * xi) * nu;
	x = *cn * *cn;
	y = *sn * *sn;
	z__ = *dn * *dn;
	w = eta * *sn;
	*b = *b + bv - w;
	*d__ = *d__ + dv + w;
	t = w / (1. - n * (y - eta * *cn * *dn));
	*j = *j + jv + uatan_(&t, &h__);
	*f = (*gel)(&u, nc, mc, b, d__, j, sn, cn, dn);
	if (v * v < *utol || fabs(*f) < *ftol) {
	    ret_val = u;
	    return ret_val;
	}
    }
    //s_wsle(&io___36);
    //do_lio(&c__9, &c__1, "(nigel) No convergence", (ftnlen)22);
    //e_wsle();
    return ret_val;
} /* nigel_ */

doublereal aigel_(doublereal *nc, doublereal *mc, doublereal *cb, doublereal *
	cd, doublereal *cj, doublereal *utol, doublereal *ftol, D_fp gel, 
	doublereal *b, doublereal *d__, doublereal *j, doublereal *sn, 
	doublereal *cn, doublereal *dn, doublereal *f)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    //integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
    //	    e_wsle(void);

    /* Local variables */
    static doublereal h__;
    static integer i__;
    static doublereal m, n, p, q, r__, t, u, v, w, x, y, z__, bh, dh, kc, bl, 
	    dl, fl, jl, bu, jh, du, fh, fu, uh, ju, xh, ul, yh, zh, yl, xi, 
	    nu, uu, yu, cnh, dnh, eta, cnl, dnl, cnu, dnu, snh, snl, snu, 
	    zeta;
    extern doublereal uatan_(doublereal *, doublereal *);

    /* Fortran I/O blocks */
    //static cilist io___85 = { 0, 6, 0, 0, 0 };


    n = 1. - *nc;
    kc = sqrt(*mc);
    m = 1. - *mc;
    h__ = n * *nc * (*mc - *nc);
    ul = 0.;
    bl = 0.;
    dl = 0.;
    jl = 0.;
    snl = 0.;
    cnl = 1.;
    dnl = 1.;
    yl = 0.;
    fl = (*gel)(&ul, nc, mc, &bl, &dl, &jl, &snl, &cnl, &dnl);
    uu = *cb + *cd;
    bu = *cb;
    du = *cd;
    ju = *cj;
    snu = 1.;
    cnu = 0.;
    dnu = kc;
    yu = 1.;
    fu = (*gel)(&uu, nc, mc, &bu, &du, &ju, &snu, &cnu, &dnu);
    uh = uu;
    bh = bu;
    dh = du;
    jh = ju;
    snh = snu;
    cnh = cnu;
    dnh = dnu;
    yh = yu;
    fh = fu;
    for (i__ = 1; i__ <= 60; ++i__) {
	uh *= .5;
	y = yh;
	v = cnh * dnh;
	p = cnh + dnh;
	q = 1. / (cnh + 1.);
	r__ = 1. / (dnh + 1.);
	xh = p * r__;
	yh = yh * q * r__;
	zh = p * q;
	w = yh * snh;
	cnh = sqrt(xh);
	snh = sqrt(yh);
	dnh = sqrt(zh);
	bh = (bh + w) * .5;
	dh = (dh - w) * .5;
	t = w / (1. - n * (y - yh * v));
	jh = (jh - uatan_(&t, &h__)) * .5;
	u = ul + uh;
	xi = cnl * cnh;
	eta = snl * snh;
	zeta = dnl * dnh;
	nu = 1. / (1. - m * yl * yh);
	*cn = (xi - eta * zeta) * nu;
	*sn = (snl * cnh * dnh + snh * cnl * dnl) * nu;
	*dn = (zeta - m * eta * xi) * nu;
	x = *cn * *cn;
	y = *sn * *sn;
	z__ = *dn * *dn;
	w = eta * *sn;
	*b = bl + bh - w;
	*d__ = dl + dh + w;
	t = w / (1. - n * (y - eta * *cn * *dn));
	*j = jl + jh + uatan_(&t, &h__);
	*f = (*gel)(&u, nc, mc, b, d__, j, sn, cn, dn);
	if (*f < 0.) {
	    ul = u;
	    bl = *b;
	    dl = *d__;
	    jl = *j;
	    snl = *sn;
	    cnl = *cn;
	    dnl = *dn;
	    yl = y;
	    fl = *f;
	} else {
	    uu = u;
	    bu = *b;
	    du = *d__;
	    ju = *j;
	    snu = *sn;
	    cnu = *cn;
	    dnu = *dn;
	    yu = y;
	    fu = *f;
	}
	if ((d__1 = uu - ul, fabs(d__1)) < *utol || fabs(*f) < *ftol) {
	    if (fabs(fu) < fabs(fl)) {
		u = uu;
		*b = bu;
		*d__ = du;
		*j = ju;
		*sn = snu;
		*cn = cnu;
		*dn = dnu;
		*f = fu;
	    } else {
		u = ul;
		*b = bl;
		*d__ = dl;
		*j = jl;
		*sn = snl;
		*cn = cnl;
		*dn = dnl;
		*f = fl;
	    }
	    ret_val = u;
	    return ret_val;
	}
    }
    //s_wsle(&io___85);
    //do_lio(&c__9, &c__1, "(aigel) No convergence", (ftnlen)22);
    //e_wsle();
    return ret_val;
} /* aigel_ */

doublereal uatan_(doublereal *t, doublereal *h__)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal), atan(doublereal), log(doublereal);

    /* Local variables */
    static doublereal r__, y, z__;

    z__ = -(*h__) * *t * *t;
    if (fabs(z__) < .001) {
	ret_val = *t * (z__ * (z__ * (z__ * (z__ * .1111111111111111 + 
		.14285714285714285) + .20000000000000001) + 
		.33333333333333331) + 1.);
    } else if (z__ < 0.) {
	r__ = sqrt(*h__);
	ret_val = atan(r__ * *t) / r__;
    } else {
	r__ = sqrt(-(*h__));
	y = r__ * *t;
	ret_val = log((y + 1.) / (1. - y)) * .5 / r__;
    }
    return ret_val;
} /* uatan_ */

/* Subroutine */ int sersdj_(doublereal *u, doublereal *n, doublereal *m, 
	doublereal *sn, doublereal *d__, doublereal *j)
{
    static doublereal d2, m2, n2, n3, n4, m3, m4, m5, m6, s1, s2, s3, s4, s5, 
	    s6, d3, d4, d5, d6, j2, j3, j4, j5, j6, u2, u3, u5, mp, mm2, mm3, 
	    mm4, m2p, m3p, m4p, m5p, m6p, mm5, m2m3, m2m4;

    n2 = *n * *n;
    n3 = n2 * *n;
    n4 = n2 * n2;
    m2 = *m * *m;
    m3 = m2 * *m;
    m4 = m2 * m2;
    m5 = m3 * m2;
    m6 = m3 * m3;
    mp = *m + 1.;
    m2p = m2 + 1.;
    m3p = m3 + 1.;
    m4p = m4 + 1.;
    m5p = m5 + 1.;
    m6p = m6 + 1.;
    mm2 = *m + m2;
    mm3 = *m + m3;
    mm4 = *m + m4;
    mm5 = *m + m5;
    m2m3 = m2 + m3;
    m2m4 = m2 + m4;
    s1 = mp * .16666666666666666;
    s2 = (m2p + *m * 14.) * .0083333333333333332;
    s3 = (m3p + mm2 * 135.) * 1.9841269841269841e-4;
    s4 = (m4p + mm3 * 1228. + m2 * 5478.) * 2.7557319223985893e-6;
    s5 = (m5p + mm4 * 11069. + m2m3 * 165826.) * 2.505210838544172e-8;
    s6 = (m6p + mm5 * 99642. + m2m4 * 4494351. + m3 * 13180268.) * 
	    1.6059043836821613e-10;
    d2 = mp * .066666666666666666;
    d3 = (m2p * 2. + *m * 13.) * .0031746031746031746;
    d4 = (m3p + mm2 * 30.) * 3.5273368606701942e-4;
    d5 = (m4p * 2. + mm3 * 251. + m2 * 876.) * 6.4133397466730804e-6;
    d6 = (m5p * 2. + mm4 * 1018. + m2m3 * 9902.) * 1.6446895031804183e-7;
    j2 = .20000000000000001;
    j3 = (mp * 30. - *n * 45.) * .0031746031746031746;
    j4 = (m2p * 63. + *m * 252. - mp * 315. * *n + n2 * 315.) * 
	    3.5273368606701942e-4;
    j5 = (m3p * 510. + mm2 * 5850. - (m2p * 6615. + *m * 21735.) * *n + mp * 
	    18900. * n2 - n3 * 14175.) * 6.4133397466730804e-6;
    j6 = (m4p * 2046. + mm3 * 59268. + m2 * 158103. - (m3p * 63360. + mm2 * 
	    497475.) * *n + (m2p * 395010. + *m * 1164240.) * n2 - mp * 
	    779625. * n3 + n4 * 467775.) * 1.6446895031804183e-7;
    u2 = *u * *u;
    u3 = u2 * *u;
    u5 = u3 * u2;
    *sn = *u * (1. - u2 * (s1 - u2 * (s2 - u2 * (s3 - u2 * (s4 - u2 * (s5 - 
	    u2 * s6))))));
    *d__ = u3 * (.33333333333333331 - u2 * (d2 - u2 * (d3 - u2 * (d4 - u2 * (
	    d5 - u2 * d6)))));
    *j = *d__ + *n * u5 * (j2 - u2 * (j3 - u2 * (j4 - u2 * (j5 - u2 * j6))));
    return 0;
} /* sersdj_ */
#if 0
doublereal bigel_(doublereal *nc, doublereal *mc, doublereal *ptol, 
	doublereal *ftol, D_fp gel, doublereal *b, doublereal *d__, 
	doublereal *j, doublereal *sn, doublereal *cn, doublereal *dn, 
	doublereal *f)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal), sqrt(doublereal);
    //integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
    //	    e_wsle(void);

    /* Local variables */
    static integer i__;
    static doublereal p, u, b0, d0, f0, j0, p0, dp, cn0, dn0, sn0;
    extern /* Subroutine */ int elbdj_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    //static cilist io___140 = { 0, 6, 0, 0, 0 };


    dp = 1.5707963267948966;
    p = dp;
    p0 = 0.;
    b0 = 0.;
    d0 = 0.;
    j0 = 0.;
    sn0 = 0.;
    cn0 = 1.;
    dn0 = 1.;
    f0 = (*gel)(&c_b15, nc, mc, &b0, &d0, &j0, &sn0, &cn0, &dn0);
    for (i__ = 1; i__ <= 60; ++i__) {
	elbdj_(&p, nc, mc, b, d__, j);
	u = *b + *d__;
	*sn = sin(p);
	*cn = cos(p);
	*dn = sqrt(*cn * *cn + *mc * *sn * *sn);
	*f = (*gel)(&u, nc, mc, b, d__, j, sn, cn, dn);
	if ((d__1 = p - p0, fabs(d__1)) < *ptol || fabs(*f) < *ftol) {
	    if (fabs(*f) > fabs(f0)) {
		p = p0;
		*b = b0;
		*d__ = d0;
		*j = j0;
		*sn = sn0;
		*cn = cn0;
		*dn = dn0;
		*f = f0;
	    }
	    ret_val = p;
	    return ret_val;
	}
	p0 = p;
	b0 = *b;
	d0 = *d__;
	j0 = *j;
	sn0 = *sn;
	cn0 = *cn;
	dn0 = *dn;
	f0 = *f;
	dp *= .5;
	if (*f < 0.) {
	    p += dp;
	} else {
	    p -= dp;
	}
    }
    //s_wsle(&io___140);
    //do_lio(&c__9, &c__1, "(bigel) No convergence", (ftnlen)22);
    //e_wsle();
    return ret_val;
} /* bigel_ */

/* Subroutine */ int elbdj_(doublereal *phi, doublereal *nc, doublereal *mc, 
	doublereal *b, doublereal *d__, doublereal *j)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sin(doublereal);

    /* Local variables */
    static doublereal f;
    extern doublereal rd_(doublereal *, doublereal *, doublereal *), rf_(
	    doublereal *, doublereal *, doublereal *), rj_(doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal sn, cn2, dn2, sn2, sn33;

    sn = sin(*phi);
    sn2 = sn * sn;
    cn2 = 1. - sn2;
    dn2 = cn2 + *mc * sn2;
    sn33 = sn2 * sn / 3.;
/* write(*,'(a10,1pe25.15,1pe25.15)') "rf:",cn2,dn2 */
    f = sn * rf_(&cn2, &dn2, &c_b20);
    *d__ = sn33 * rd_(&cn2, &dn2, &c_b20);
    d__1 = cn2 + *nc * sn2;
    *j = sn33 * rj_(&cn2, &dn2, &c_b20, &d__1);
    *b = f - *d__;
    return 0;
} /* elbdj_ */
#endif
