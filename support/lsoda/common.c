
int      g_nyh = 0, g_lenyh = 0;

/* newly added static variables */

int      imxer;
int      mord[3] = {0, 12, 5};
double   sqrteta, *yp1, *yp2;
double   sm1[13] = {0., 0.5, 0.575, 0.55, 0.45, 0.35, 0.25, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025};

/* static variables for lsoda() */

double   ccmax, el0, h, hu, rc, tn;
int      illin = 0, init = 0, nhnil, ntrep = 0, nslast, nyh, ierpj, iersl,
                jcur, l, meth, miter, maxord, maxcor, msbp, mxncf, n, nq, nst,
                nfe, nje, nqu;
double   tsw, pdnorm;
int      jtyp, mused;

/* no static variable for prja(), solsy() */
/* static variables for stoda() */

double   conit, crate, el[14], elco[13][14], hold, rmax, tesco[13][4];
int      ialth, ipup, lmax, nslp;
double   pdest, pdlast, ratio, cm1[13], cm2[6];
int      icount, irflag;

/* static variables for various vectors and the Jacobian. */

double **yh, **wm, *ewt, *savf, *acor;
int     *ipvt;

