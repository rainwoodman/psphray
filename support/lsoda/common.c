#include "lsoda.h"
#include "lsoda_internal.h"

/* newly added static variables */

int      imxer;
int      mord[3] = {0, 12, 5};
double   sm1[13] = {0., 0.5, 0.575, 0.55, 0.45, 0.35, 0.25, 0.2, 0.15, 0.1, 0.075, 0.05, 0.025};

/* static variables for lsoda() */

double   h, hu, rc, tn;
int      illin = 0, init = 0, nhnil, ntrep = 0, nslast,
                jcur, meth, mused, nq, nst,
                nfe, nje, nqu, miter;
double   tsw, pdnorm;

/* no static variable for prja(), solsy() */
/* static variables for stoda() */

double   crate, el[14], elco[13][14], hold, rmax, tesco[13][4];
int      ialth, ipup, nslp;
double   pdest, pdlast, cm1[13], cm2[6];
int      icount, irflag;

/* static variables for various vectors and the Jacobian. */

void * memory;

double **yh, **wm, *ewt, *savf, *acor;
int     *ipvt;