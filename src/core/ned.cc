/* 
*/

#ifdef __cplusplus
extern "C" {
#endif


#include<math.h>


/* Table of constant values */

double max(double d__1, double d__2)
{
  if (d__1>d__2) return(d__1); return(d__2);
};

double min(double d__1, double d__2)
{
  if (d__1>d__2) return(d__2); return(d__1);
};

long min_ii(long d__1, long d__2)
{
  if (d__1>d__2) return(d__2); return(d__1);
};

long max_ii(long d__1, long d__2)
{
  if (d__1>d__2) return(d__1); return(d__2);
};


static long c_n1 = -1;

/* 	ARTURO QUIRANTES SIERRA */
/* 	Department of Applied Physics, Faculty of Sciences */
/* 	University of Granada, 18071 Granada (SPAIN) */
/* 	http://www.ugr.es/local/aquiran/codigos.htm */
/* 	aquiran@ugr.es */

/* 	Last update: 20 May 2.003 */

/* 	Subroutine NED */
/* 	to calculate Clebsch-Gordan coefficients */

/* 	You need to add a "NED(AJ,BJ,CJ,AM,BM,CM,CG)" in your main routine */
/* 	Input: */
/* 		AJ,BJ,CJ,AM,BM,CM (the usual Clebsch-Gordan indices) */
/* 	Output: */
/* 		CG=C-G(AJ,BJ,CJ,AM,BM,CM) */

/* Subroutine */ 
long ned_(double *aj, double *bj, double *cj, 
	double *am, double *bm, double *cm, double *cg)
{
    /* System generated locals */
    long i__1, i__2;
    double d__1, d__2;

    /* Builtin functions */
//    double sqrt(double);
//    long i_dnnt(double *), pow_ii(long *, long *);

    /* Local variables */
    static double d__;
    static long i__, k;
    static double p, q[10000]	/* was [100][100] */, x;
    static long i2, k0, k1, ja, ma, jb, mb, jc, mc, la, lb, lc, ld;
    static double fn;
    static long ip, lt, zz, ja2, jb2, jc2;

/* Computing MAX */
    d__1 = *aj * 2 + 1, d__2 = *bj * 2 + 1, d__1 = max(d__1,d__2), d__2 = *cj 
	    * 2 + 1, d__1 = max(d__1,d__2), d__2 = *aj + *bj + *cj, d__1 = 
	    max(d__1,d__2), d__2 = *aj + *am, d__1 = max(d__1,d__2), d__2 = *
	    bj + *bm, d__1 = max(d__1,d__2), d__2 = *cj + *cm;
    zz = (long) (max(d__1,d__2) + 2);
    i__1 = zz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ - 1] = 1.;
	q[i__ + i__ * 100 - 101] = 1.;
/* L2: */
    }
    i__1 = zz - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = i__;
	for (k = 2; k <= i__2; ++k) {
	    q[i__ + 1 + k * 100 - 101] = q[i__ + (k - 1) * 100 - 101] + q[i__ 
		    + k * 100 - 101];
/* L3: */
	}
    }
    *cg = 0.;
    ja = (long) (*aj + *am + 1.01);
    ma = (long) (*aj - *am + 1.01);
    jb = (long) (*bj + *bm + 1.01);
    mb = (long) (*bj - *bm + 1.01);
    jc = (long) (*cj + *cm + 1.01);
    mc = (long) (*cj - *cm + 1.01);
    la = (long) (*bj + *cj - *aj + 1.01);
    lb = (long) (*cj + *aj - *bj + 1.01);
    lc = (long) (*aj + *bj - *cj + 1.01);
    lt = (long) (*aj + *bj + *cj + 1.01);
    d__ = (d__1 = *am + *bm - *cm, fabs(d__1)) - .01;
    if (d__ <= 0.) {
	goto L10;
    } else {
	goto L20;
    }
L10:
/* Computing MIN */
    i__2 = min_ii(ja,jb), i__2 = min_ii(i__2,jc), i__2 = min_ii(i__2,ma),
    i__2 = min_ii(i__2,mb), i__2 = min_ii(i__2,mc), i__2 = min_ii(i__2,la),
    i__2 = min_ii(i__2,lb);
    ld   = min_ii(i__2,lc);
    if (ld <= 0) {
	goto L20;
    } else {
	goto L30;
    }
L30:
    ja2 = (long) (*aj + *aj + *am + *am);
    jb2 = (long) (*bj + *bj + *bm + *bm);
    jc2 = (long) (*cj + *cj - *cm - *cm);
    i2 = ja2 + jb2 + jc2 - (ja2 / 2 << 1) - (jb2 / 2 << 1) - (jc2 / 2 << 1);
    if (i2 != 0) {
	goto L20;
    } else {
	goto L40;
    }
L40:
    fn = q[ja + ma - 1 + lc * 100 - 101] / q[lt + (jc + mc - 1) * 100 - 101];
    fn = fn * q[jb + mb - 1 + lc * 100 - 101] / q[lt + 100];
    fn /= q[ja + ma - 1 + ja * 100 - 101];
    fn /= q[jb + mb - 1 + jb * 100 - 101];
    fn /= q[jc + mc - 1 + jc * 100 - 101];
/* Computing MAX */
    i__2 = 0, i__1 = lc - ja, i__2 = max_ii(i__2,i__1), i__1 = lc - mb;
    k0 = max_ii(i__2,i__1) + 1;
/* Computing MIN */
    i__2 = min_ii(lc,ma);
    k1 = min_ii(i__2,jb);
    x = 0.;
    i__2 = k1;
    for (k = k0; k <= i__2; ++k) {
	x = -x - q[lc + k * 100 - 101] * q[lb + (ma - k + 1) * 100 - 101] * q[
		la + (jb - k + 1) * 100 - 101];
/* L50: */
    }
    ip = k1 + lb + jc;
    p = (double) (1 - (ip - (ip / 2 << 1) << 1));
    *cg = p * x * sqrt(fn);
/* 	What we?ve calculated is a Wigner 3-j coefficient */
/* 	Next, we?ll turn into a Clebsch-Gordan coefficient */
    d__1 = *aj - *bj - *cm;
    i__2 = (long)(d__1);
    *cg = *cg * sqrt(*cj * 2 + 1) * pow(-1.0, i__2);
L20:
    return 0;
} /* ned_ */

#ifdef __cplusplus
	}
#endif
