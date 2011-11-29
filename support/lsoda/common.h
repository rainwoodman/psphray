
#define max( a , b )  ( (a) > (b) ? (a) : (b) )
#define min( a , b )  ( (a) < (b) ? (a) : (b) )

#define ETA 2.2204460492503131e-16
#define SQRTETA 1.4901161193847656e-08
#define CCMAX  0.3
#define MAXCOR 3
#define MSBP 20
#define MXNCF 10
#define RATIO 5.0

extern const int      mord[3];
extern const double   sm1[13];
/* newly added static variables */

struct common_t {
	int      imxer;

	/* static variables for lsoda() */

	double   h, hu, rc, tn;
	int      illin, init, nhnil, ntrep, nslast,
					jcur, meth, mused, nq, nst,
					nfe, nje, nqu, miter;
	double   tsw, pdnorm;

	/* no static variable for prja(), solsy() */
	/* static variables for stoda() */

	double   crate, el[14];
#ifdef CFODE_STATIC
	double (*elco)[14], (*tesco)[4];
#else
	double elco[13][14], tesco[13][4];
#endif
	double hold, rmax;
	int      ialth, ipup, nslp;
	double   pdest, pdlast, cm1[13], cm2[6];
	int      icount, irflag;

	/* static variables for various vectors and the Jacobian. */

	void * memory;
	double **yh, **wm, *ewt, *savf, *acor;
	int     *ipvt;
};
#define _C(x) (common->x)
