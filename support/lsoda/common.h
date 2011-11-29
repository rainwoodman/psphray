
#define max( a , b )  ( (a) > (b) ? (a) : (b) )
#define min( a , b )  ( (a) < (b) ? (a) : (b) )

#define ETA 2.2204460492503131e-16
#define SQRTETA 1.4901161193847656e-08
#define CCMAX  0.3
#define MAXCOR 3
#define MSBP 20
#define MXNCF 10
#define RATIO 5.0

/* newly added static variables */

extern int      imxer;
extern int      mord[3];
extern double   sm1[13];

/* static variables for lsoda() */

extern double   h, hu, rc, tn;
extern int      illin, init, nhnil, ntrep, nslast,
                jcur, meth, mused, nq, nst,
                nfe, nje, nqu, miter;
extern double   tsw, pdnorm;

/* no static variable for prja(), solsy() */
/* static variables for stoda() */

extern double   crate, el[14], elco[13][14], hold, rmax, tesco[13][4];
extern int      ialth, ipup, nslp;
extern double   pdest, pdlast, cm1[13], cm2[6];
extern int      icount, irflag;

/* static variables for various vectors and the Jacobian. */

extern void * memory;
extern struct vec_t vec;
