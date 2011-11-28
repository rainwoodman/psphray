/*
  This is a C version of the LSODA library. I acquired the original
  source code from this web page:

    http://www.ccl.net/cca/software/SOURCES/C/kinetics2/index.shtml

  I merged several C files into one and added a simpler interface. I
  also made the array start from zero in functions called by lsoda(),
  and fixed two minor bugs: a) small memory leak in freevectors(); and
  b) misuse of lsoda() in the example.

  The original source code came with no license or copyright
  information. I now release this file under the MIT/X11 license. All
  authors' notes are kept in this file.

  - Heng Li <lh3lh3@gmail.com>
 */

/* The MIT License

   Copyright (c) 2009 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdio.h>
#include <math.h>
#include "lsoda.h"
#include "common.h"

/***********
 * lsoda.c *
 ***********/

/*
From tam@dragonfly.wri.com Wed Apr 24 01:35:52 1991
Return-Path: <tam>
Date: Wed, 24 Apr 91 03:35:24 CDT
From: tam@dragonfly.wri.com
To: whitbeck@wheeler.wrc.unr.edu
Subject: lsoda.c
Cc: augenbau@sparc0.brc.uconn.edu


I'm told by Steve Nichols at Georgia Tech that you are interested in
a stiff integrator.  Here's a translation of the fortran code LSODA.

Please note
that there is no comment.  The interface is the same as the FORTRAN
code and I believe the documentation in LSODA will suffice.
As usual, a free software comes with no guarantee.

Hon Wah Tam
Wolfram Research, Inc.
tam@wri.com
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lsoda_internal.h"
#include "blas.h"

static void     terminate(int *istate);
static void     terminate2(double *y, double *t);
static void     successreturn(double *y, double *t, int itask, int ihit, double tcrit, int *istate);
/* Terminate lsoda due to illegal input. */
static void terminate(int *istate)
{
	if (illin == 5) {
		fprintf(stderr, "[lsoda] repeated occurrence of illegal input. run aborted.. apparent infinite loop\n");
	} else {
		illin++;
		*istate = -3;
	}
}


/* Terminate lsoda due to various error conditions. */
static void terminate2(double *y, double *t)
{
	int             i;
	double ** yh = vec.yh;
	for (i = 1; i <= n; i++)
		y[i] = yh[1][i];
	*t = tn;
	illin = 0;
	return;

}

/*
   The following block handles all successful returns from lsoda.
   If itask != 1, y is loaded from yh and t is set accordingly.
   *Istate is set to 2, the illegal input counter is zeroed, and the
   optional outputs are loaded into the work arrays before returning.
*/

static void successreturn(double *y, double *t, int itask, int ihit, double tcrit, int *istate)
{
	int             i;
	double ** yh = vec.yh;
	for (i = 1; i <= n; i++)
		y[i] = yh[1][i];
	*t = tn;
	if (itask == 4 || itask == 5)
		if (ihit)
			*t = tcrit;
	*istate = 2;
	illin = 0;
}

static int check_opt(struct lsoda_opt_t * opt, int *istate, int n) {
	const int mxstp0 = 500, mxhnl0 = 10;

	if (*istate == 1) {
		opt->h0 = 0.;
		opt->mxordn = mord[1];
		opt->mxords = mord[2];
	}
	if (opt->ixpr < 0 || opt->ixpr > 1) {
		fprintf(stderr, "[lsoda] ixpr = %d is illegal\n", opt->ixpr);
		terminate(istate);
		return 0;
	}
	if (opt->mxstep < 0) {
		fprintf(stderr, "[lsoda] mxstep < 0\n");
		terminate(istate);
		return 0;
	}
	if (opt->mxstep == 0) opt->mxstep = mxstp0;
	if (opt->mxhnil < 0) {
		fprintf(stderr, "[lsoda] mxhnil < 0\n");
		terminate(istate);
		return 0;
	}
	if (*istate == 1) {
		if (opt->mxordn < 0) {
			fprintf(stderr, "[lsoda] mxordn = %d is less than 0\n", opt->mxordn);
			terminate(istate);
			return 0;
		}
		if (opt->mxordn == 0) opt->mxordn = 100;
		opt->mxordn = min(opt->mxordn, mord[1]);
		if (opt->mxords < 0) {
			fprintf(stderr, "[lsoda] mxords = %d is less than 0\n", opt->mxords);
			terminate(istate);
			return 0;
		}
		if (opt->mxords == 0) opt->mxords = 100;
		opt->mxords = min(opt->mxords, mord[2]);
	}	/* end if ( *istate == 1 )  */
	if (opt->hmax < 0.) {
		fprintf(stderr, "[lsoda] hmax < 0.\n");
		terminate(istate);
		return 0;
	}
	opt->hmxi = 0.;
	if (opt->hmax > 0)
		opt->hmxi = 1. / opt->hmax;
	if (opt->hmin < 0.) {
		fprintf(stderr, "[lsoda] hmin < 0.\n");
		terminate(istate);
		return 0;
	}
	return 1;
}
static int alloc_mem(struct lsoda_opt_t * opt, int n) {
	int nyh = n;
	int lenyh = 1 + max(opt->mxordn, opt->mxords);
	long offset = 0;
	int i;
	long yhoff = offset;
	/* yh */
	offset += (1 + lenyh) * sizeof(double *);
	long yh0off = offset;
	for(i = 0; i <= lenyh; i++) {
		offset += (1 + nyh) * sizeof(double);
	}

	long wmoff = offset;
	long wm0off = offset;
	/* wm */
	offset += (1 + nyh) * sizeof(double *);
	for(i = 0; i <= nyh; i++) {
		offset += (1 + nyh) * sizeof(double);
	}
	
	/* ewt */
	long ewtoff =  offset;
	offset += (1 + nyh) * sizeof(double);

	/* savf */
	long savfoff = offset;
	offset += (1 + nyh) * sizeof(double);

	/* acor */
	long acoroff = offset;
	offset += (1 + nyh) * sizeof(double);

	/* ipvt */
	long ipvtoff = offset;
	offset += (1 + nyh) * sizeof(int);

	memory = malloc(offset);

	vec.yh = memory + yhoff;
	vec.wm =  memory + wmoff;
	vec.ewt = memory + ewtoff;
	vec.savf =memory + savfoff;
	vec.acor =memory + acoroff;
	vec.ipvt =memory + ipvtoff;

	for(i = 0; i <= lenyh; i++) {
		vec.yh[i] = memory + yh0off + i * (1 + nyh) * sizeof(double);
	}

	for(i = 0; i <= nyh; i++) {
		vec.wm[i] = memory + wm0off + i * (1 + nyh) * sizeof(double);
	}

	return memory != NULL;
}
/*
c-----------------------------------------------------------------------
c this is the march 30, 1987 version of
c lsoda.. livermore solver for ordinary differential equations, with
c         automatic method switching for stiff and nonstiff problems.
c
c this version is in double precision.
c
c lsoda solves the initial value problem for stiff or nonstiff
c systems of first order ode-s,
c     dy/dt = f(t,y) ,  or, in component form,
c     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
c
c this a variant version of the lsode package.
c it switches automatically between stiff and nonstiff methods.
c this means that the user does not have to determine whether the
c problem is stiff or not, and the solver will automatically choose the
c appropriate method.  it always starts with the nonstiff method.
c
c authors..
c                linda r. petzold  and  alan c. hindmarsh,
c                computing and mathematics research division, l-316
c                lawrence livermore national laboratory
c                livermore, ca 94550.
c
c references..
c 1.  alan c. hindmarsh,  odepack, a systematized collection of ode
c     solvers, in scientific computing, r. s. stepleman et al. (eds.),
c     north-holland, amsterdam, 1983, pp. 55-64.
c 2.  linda r. petzold, automatic selection of methods for solving
c     stiff and nonstiff systems of ordinary differential equations,
c     siam j. sci. stat. comput. 4 (1983), pp. 136-148.
c-----------------------------------------------------------------------
c summary of usage.
c
c communication between the user and the lsoda package, for normal
c situations, is summarized here.  this summary describes only a subset
c of the full set of options available.  see the full description for
c details, including alternative treatment of the jacobian matrix,
c optional inputs and outputs, nonstandard options, and
c instructions for special situations.  see also the example
c problem (with program and output) following this summary.
c
c a. first provide a subroutine of the form..
c               subroutine f (neq, t, y, ydot)
c               dimension y(neq), ydot(neq)
c which supplies the vector function f by loading ydot(i) with f(i).
c
c b. write a main program which calls subroutine lsoda once for
c each point at which answers are desired.  this should also provide
c for possible use of logical unit 6 for output of error messages
c by lsoda.  on the first call to lsoda, supply arguments as follows..
c f      = name of subroutine for right-hand side vector f.
c          this name must be declared external in calling program.
c neq    = number of first order ode-s.
c y      = array of initial values, of length neq.
c t      = the initial value of the independent variable.
c tout   = first point where output is desired (.ne. t).
c itol   = 1 or 2 according as atol (below) is a scalar or array.
c rtol   = relative tolerance parameter (scalar).
c atol   = absolute tolerance parameter (scalar or array).
c          the estimated local error in y(i) will be controlled so as
c          to be less than
c             ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or
c             ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2.
c          thus the local error test passes if, in each component,
c          either the absolute error is less than atol (or atol(i)),
c          or the relative error is less than rtol.
c          use rtol = 0.0 for pure absolute error control, and
c          use atol = 0.0 (or atol(i) = 0.0) for pure relative error
c          control.  caution.. actual (global) errors may exceed these
c          local tolerances, so choose them conservatively.
c itask  = 1 for normal computation of output values of y at t = tout.
c istate = integer flag (input and output).  set istate = 1.
c iopt   = 0 to indicate no optional inputs used.
c rwork  = real work array of length at least..
c             22 + neq * max(16, neq + 9).
c          see also paragraph e below.
c lrw    = declared length of rwork (in user-s dimension).
c iwork  = integer work array of length at least  20 + neq.
c liw    = declared length of iwork (in user-s dimension).
c jac    = name of subroutine for jacobian matrix.
c          use a dummy name.  see also paragraph e below.
c jt     = jacobian type indicator.  set jt = 2.
c          see also paragraph e below.
c note that the main program must declare arrays y, rwork, iwork,
c and possibly atol.
c
c c. the output from the first call (or any call) is..
c      y = array of computed values of y(t) vector.
c      t = corresponding value of independent variable (normally tout).
c istate = 2  if lsoda was successful, negative otherwise.
c          -1 means excess work done on this call (perhaps wrong jt).
c          -2 means excess accuracy requested (tolerances too small).
c          -3 means illegal input detected (see printed message).
c          -4 means repeated error test failures (check all inputs).
c          -5 means repeated convergence failures (perhaps bad jacobian
c             supplied or wrong choice of jt or tolerances).
c          -6 means error weight became zero during problem. (solution
c             component i vanished, and atol or atol(i) = 0.)
c          -7 means work space insufficient to finish (see messages).
c
c d. to continue the integration after a successful return, simply
c reset tout and call lsoda again.  no other parameters need be reset.
c
c e. note.. if and when lsoda regards the problem as stiff, and
c switches methods accordingly, it must make use of the neq by neq
c jacobian matrix, j = df/dy.  for the sake of simplicity, the
c inputs to lsoda recommended in paragraph b above cause lsoda to
c treat j as a full matrix, and to approximate it internally by
c difference quotients.  alternatively, j can be treated as a band
c matrix (with great potential reduction in the size of the rwork
c array).  also, in either the full or banded case, the user can supply
c j in closed form, with a routine whose name is passed as the jac
c argument.  these alternatives are described in the paragraphs on
c rwork, jac, and jt in the full description of the call sequence below.
c
c-----------------------------------------------------------------------
*/

void lsoda(_lsoda_f f, int neq, double *y, double *t, double tout, int itol, double *rtol, double *atol,
		int itask, int *istate, struct lsoda_opt_t * opt, void *_data) {

		int kflag;
		int jstart;
		struct lsoda_opt_t def_opt = {0};

		int             i, iflag, lenyh, ihit;
		double          atoli, ayi, big, h0, hmx, rh, rtoli, tcrit, tdist, tnext, tol,
						tolsf, tp, size, sum, w0;

		/*
		   Block a.
		   This code block is executed on every call.
		   It tests *istate and itask for legality and branches appropriately.
		   If *istate > 1 but the flag init shows that initialization has not
		   yet been done, an error return occurs.
		   If *istate = 1 and tout = t, return immediately.
		 */

		if (*istate < 1 || *istate > 3) {
			fprintf(stderr, "[lsoda] illegal istate = %d\n", *istate);
			terminate(istate);
			return;
		}
		if (itask < 1 || itask > 5) {
			fprintf(stderr, "[lsoda] illegal itask = %d\n", itask);
			terminate(istate);
			return;
		}
		if (init == 0 && (*istate == 2 || *istate == 3)) {
			fprintf(stderr, "[lsoda] istate > 1 but lsoda not initialized\n");
			terminate(istate);
			return;
		}
		if (*istate == 1) {
			init = 0;
			if (tout == *t) {
				ntrep++;
				if (ntrep < 5) return;
				fprintf(stderr, "[lsoda] repeated calls with istate = 1 and tout = t. run aborted.. apparent infinite loop\n");
				return;
			}
		}
		/*
		   Block b.
		   The next code block is executed for the initial call ( *istate = 1 ),
		   or for a continuation call with parameter changes ( *istate = 3 ).
		   It contains checking of all inputs and various initializations.

		   First check legality of the non-optional inputs neq, itol, iopt,
		   jt, ml, and mu.
		 */

		if (*istate == 1 || *istate == 3) {
			ntrep = 0;
			if (neq <= 0) {
				fprintf(stderr, "[lsoda] neq = %d is less than 1\n", neq);
				terminate(istate);
				return;
			}
			if (*istate == 3 && neq > n) {
				fprintf(stderr, "[lsoda] istate = 3 and neq increased\n");
				terminate(istate);
				return;
			}
			n = neq;
			if (itol < 1 || itol > 4) {
				fprintf(stderr, "[lsoda] itol = %d illegal\n", itol);
				terminate(istate);
				return;
			}
			/* Next process and check the optional inpus.   */

			/* Default options.   */

			if(opt == NULL) {
				fprintf(stderr, "[lsoda] need an opt struct\n");
				terminate(istate);
				return;
			}
			if(!check_opt(opt, istate, n)) {
				return;
			}
			h0 = opt->h0;
			if(*istate == 1) {
				if ((tout - *t) * h0 < 0.) {
					fprintf(stderr, "[lsoda] tout = %g behind t = %g. integration direction is given by %g\n",
							tout, *t, h0);
					terminate(istate);
					return;
				}
			}
		}			/* end if ( *istate == 1 || *istate == 3 )   */

		/*
		   If *istate = 1, meth is initialized to 1.

		   Also allocate memory for yh, wm, ewt, savf, acor, ipvt.
		 */
		if (*istate == 1) {
			/*
			   If memory were not freed, *istate = 3 need not reallocate memory.
			   Hence this section is not executed by *istate = 3.
			 */
			meth = 1;
			if(!alloc_mem(opt, n)) {
				printf("lsoda -- insufficient memory for your problem\n");
				terminate(istate);
				return;
			}

		}
		/*
		   Check rtol and atol for legality.
		 */
		if (*istate == 1 || *istate == 3) {
			rtoli = rtol[1];
			atoli = atol[1];
			for (i = 1; i <= n; i++) {
				if (itol >= 3)
					rtoli = rtol[i];
				if (itol == 2 || itol == 4)
					atoli = atol[i];
				if (rtoli < 0.) {
					fprintf(stderr, "[lsoda] rtol = %g is less than 0.\n", rtoli);
					terminate(istate);
					return;
				}
				if (atoli < 0.) {
					fprintf(stderr, "[lsoda] atol = %g is less than 0.\n", atoli);
					terminate(istate);
					return;
				}
			}		/* end for   */
		}			/* end if ( *istate == 1 || *istate == 3 )   */
		/*
		   If *istate = 3, set flag to signal parameter changes to stoda.
		 */
		if (*istate == 3) {
			jstart = -1;
		}
		/*
		   Block c.
		   The next block is for the initial call only ( *istate = 1 ).
		   It contains all remaining initializations, the initial call to f,
		   and the calculation of the initial step size.
		   The error weights in ewt are inverted after being loaded.
		 */
		if (*istate == 1) {
			double ** yh = vec.yh;
			double * ewt = vec.ewt;
			tn = *t;
			tsw = *t;
			maxord = opt->mxordn;
			if (itask == 4 || itask == 5) {
				tcrit = opt->tcrit;
				if ((tcrit - tout) * (tout - *t) < 0.) {
					fprintf(stderr, "[lsoda] itask = 4 or 5 and tcrit behind tout\n");
					terminate(istate);
					return;
				}
				if (h0 != 0. && (*t + h0 - tcrit) * h0 > 0.)
					h0 = tcrit - *t;
			}
			jstart = 0;
			nhnil = 0;
			nst = 0;
			nje = 0;
			nslast = 0;
			hu = 0.;
			nqu = 0;
			mused = 0;
			miter = 0;
			ccmax = 0.3;
			maxcor = 3;
			msbp = 20;
			mxncf = 10;
			/*
			   Initial call to f.
			 */
			(*f) (*t, y + 1, yh[2] + 1, _data);
			nfe = 1;
			/*
			   Load the initial value vector in yh.
			 */
			for (i = 1; i <= n; i++)
				yh[1][i] = y[i];

			/*
			   Load and invert the ewt array.  
			 */
			if(!ewset(vec.ewt, itol, rtol, atol, y, n)) {
				terminate2(y, t);
				return;
			}

			/*		( h is temporarily set to 1. ) */
			nq = 1; 
			h = 1.;

			/*
			   The coding below computes the step size, h0, to be attempted on the
			   first step, unless the user has supplied a value for this.
			   First check that tout - *t differs significantly from zero.
			   A scalar tolerance quantity tol is computed, as max(rtol[i])
			   if this is positive, or max(atol[i]/fabs(y[i])) otherwise, adjusted
			   so as to be between 100*ETA and 0.001.
			   Then the computed value h0 is given by

			   h0^(-2) = 1. / ( tol * w0^2 ) + tol * ( norm(f) )^2

			   where   w0     = max( fabs(*t), fabs(tout) ),
			   f      = the initial value of the vector f(t,y), and
			   norm() = the weighted vector norm used throughout, given by
			   the vmnorm function routine, and weighted by the
			   tolerances initially loaded into the ewt array.

			   The sign of h0 is inferred from the initial values of tout and *t.
			   fabs(h0) is made < fabs(tout-*t) in any case.
			 */
			if (h0 == 0.) {
				tdist = fabs(tout - *t);
				w0 = fmax(fabs(*t), fabs(tout));
				if (tdist < 2. * ETA * w0) {
					fprintf(stderr, "[lsoda] tout too close to t to start integration\n ");
					terminate(istate);
					return;
				}
				tol = rtol[1];
				if (itol > 2) {
					for (i = 2; i <= n; i++)
						tol = fmax(tol, rtol[i]);
				}
				if (tol <= 0.) {
					atoli = atol[1];
					for (i = 1; i <= n; i++) {
						if (itol == 2 || itol == 4)
							atoli = atol[i];
						ayi = fabs(y[i]);
						if (ayi != 0.)
							tol = fmax(tol, atoli / ayi);
					}
				}
				tol = fmax(tol, 100. * ETA);
				tol = fmin(tol, 0.001);
				sum = vmnorm(n, yh[2], ewt);
				sum = 1. / (tol * w0 * w0) + tol * sum * sum;
				h0 = 1. / sqrt(sum);
				h0 = fmin(h0, tdist);
				h0 = h0 * ((tout - *t >= 0.) ? 1. : -1.);
			}		/* end if ( h0 == 0. )   */
			/*
			   Adjust h0 if necessary to meet hmax bound.
			 */
			rh = fabs(h0) * opt->hmxi;
			if (rh > 1.)
				h0 /= rh;
			/*
			   Load h with h0 and scale yh[2] by h0.
			 */
			h = h0;
			for (i = 1; i <= n; i++)
				yh[2][i] *= h0;
		}			/* if ( *istate == 1 )   */
		/*
		   Block d.
		   The next code block is for continuation calls only ( *istate = 2 or 3 )
		   and is to check stop conditions before taking a step.
		 */
		if (*istate == 2 || *istate == 3) {
			jstart = 1;
			nslast = nst;
			switch (itask) {
				case 1:
					if ((tn - tout) * h >= 0.) {
						intdy(tout, 0, y, &iflag);
						if (iflag != 0) {
							fprintf(stderr, "[lsoda] trouble from intdy, itask = %d, tout = %g\n", itask, tout);
							terminate(istate);
							return;
						}
						*t = tout;
						*istate = 2;
						illin = 0;
						return;
					}
					break;
				case 2:
					break;
				case 3:
					tp = tn - hu * (1. + 100. * ETA);
					if ((tp - tout) * h > 0.) {
						fprintf(stderr, "[lsoda] itask = %d and tout behind tcur - hu\n", itask);
						terminate(istate);
						return;
					}
					if ((tn - tout) * h < 0.) break;
					successreturn(y, t, itask, ihit, tcrit, istate);
					return;
				case 4:
					tcrit = opt->tcrit;
					if ((tn - tcrit) * h > 0.) {
						fprintf(stderr, "[lsoda] itask = 4 or 5 and tcrit behind tcur\n");
						terminate(istate);
						return;
					}
					if ((tcrit - tout) * h < 0.) {
						fprintf(stderr, "[lsoda] itask = 4 or 5 and tcrit behind tout\n");
						terminate(istate);
						return;
					}
					if ((tn - tout) * h >= 0.) {
						intdy(tout, 0, y, &iflag);
						if (iflag != 0) {
							fprintf(stderr, "[lsoda] trouble from intdy, itask = %d, tout = %g\n", itask, tout);
							terminate(istate);
							return;
						}
						*t = tout;
						*istate = 2;
						illin = 0;
						return;
					}
				case 5:
					if (itask == 5) {
						tcrit = opt->tcrit;
						if ((tn - tcrit) * h > 0.) {
							fprintf(stderr, "[lsoda] itask = 4 or 5 and tcrit behind tcur\n");
							terminate(istate);
							return;
						}
					}
					hmx = fabs(tn) + fabs(h);
					ihit = fabs(tn - tcrit) <= (100. * ETA * hmx);
					if (ihit) {
						*t = tcrit;
						successreturn(y, t, itask, ihit, tcrit, istate);
						return;
					}
					tnext = tn + h * (1. + 4. * ETA);
					if ((tnext - tcrit) * h <= 0.)
						break;
					h = (tcrit - tn) * (1. - 4. * ETA);
					if (*istate == 2)
						jstart = -2;
					break;
			}		/* end switch   */
		}			/* end if ( *istate == 2 || *istate == 3 )   */
		/*
		   Block e.
		   The next block is normally executed for all calls and contains
		   the call to the one-step core integrator stoda.

		   This is a looping point for the integration steps.

		   First check for too many steps being taken, update ewt ( if not at
		   start of problem).  Check for too much accuracy being requested, and
		   check for h below the roundoff level in *t.
		 */
		while (1) {
			if (*istate != 1 || nst != 0) {
				if ((nst - nslast) >= opt->mxstep) {
					fprintf(stderr, "[lsoda] %d steps taken before reaching tout\n", opt->mxstep);
					*istate = -1;
					terminate2(y, t);
					return;
				}
				if(!ewset(vec.ewt, itol, rtol, atol, vec.yh[1], n)) {
					terminate2(y, t);
					*istate = -6;
					return;
				}
			}
			tolsf = ETA * vmnorm(n, vec.yh[1], vec.ewt);
			if (tolsf > 0.01) {
				tolsf = tolsf * 200.;
				if (nst == 0) {
					fprintf(stderr, "lsoda -- at start of problem, too much accuracy\n");
					fprintf(stderr, "         requested for precision of machine,\n");
					fprintf(stderr, "         suggested scaling factor = %g\n", tolsf);
					terminate(istate);
					return;
				}
				fprintf(stderr, "lsoda -- at t = %g, too much accuracy requested\n", *t);
				fprintf(stderr, "         for precision of machine, suggested\n");
				fprintf(stderr, "         scaling factor = %g\n", tolsf);
				*istate = -2;
				terminate2(y, t);
				return;
			}
			if ((tn + h) == tn) {
				nhnil++;
				if (nhnil <= opt->mxhnil) {
					fprintf(stderr, "lsoda -- warning..internal t = %g and h = %g are\n", tn, h);
					fprintf(stderr, "         such that in the machine, t + h = t on the next step\n");
					fprintf(stderr, "         solver will continue anyway.\n");
					if (nhnil == opt->mxhnil) {
						fprintf(stderr, "lsoda -- above warning has been issued %d times,\n", nhnil);
						fprintf(stderr, "         it will not be issued again for this problem\n");
					}
				}
			}
			/*
			   Call stoda
			 */
			kflag = stoda(neq, y, f, _data, jstart, opt->hmxi, opt->hmin, opt->mxords, opt->mxordn);
			/*
			   printf( "meth= %d,   order= %d,   nfe= %d,   nje= %d\n",
			   meth, nq, nfe, nje );
			   printf( "t= %20.15e,   h= %20.15e,   nst=%d\n", tn, h, nst );
			   printf( "y= %20.15e,   %20.15e,   %20.15e\n\n\n",
			   yh[1][1], yh[1][2], yh[1][3] );
			 */

			if (kflag == 0) {
				/*
				   Block f.
				   The following block handles the case of a successful return from the
				   core integrator ( kflag = 0 ).
				   If a method switch was just made, record tsw, reset maxord,
				   set jstart to -1 to signal stoda to complete the switch,
				   and do extra printing of data if ixpr = 1.
				   Then, in any case, check for stop conditions.
				 */
				jstart = 1;
				init = 1;
				if (meth != mused) {
					tsw = tn;
					maxord = opt->mxordn;
					if (meth == 2)
						maxord = opt->mxords;
					jstart = -1;
					if (opt->ixpr) {
						if (meth == 2)
							fprintf(stderr, "[lsoda] a switch to the stiff method has occurred ");
						if (meth == 1)
							fprintf(stderr, "[lsoda] a switch to the nonstiff method has occurred");
						fprintf(stderr, "at t = %g, tentative step size h = %g, step nst = %d\n", tn, h, nst);
					}
				}	/* end if ( meth != mused )   */
				/*
				   itask = 1.
				   If tout has been reached, interpolate.
				 */
				if (itask == 1) {
					if ((tn - tout) * h < 0.)
						continue;
					intdy(tout, 0, y, &iflag);
					*t = tout;
					*istate = 2;
					illin = 0;
					return;
				}
				/*
				   itask = 2.
				 */
				if (itask == 2) {
					successreturn(y, t, itask, ihit, tcrit, istate);
					return;
				}
				/*
				   itask = 3.
				   Jump to exit if tout was reached.
				 */
				if (itask == 3) {
					if ((tn - tout) * h >= 0.) {
						successreturn(y, t, itask, ihit, tcrit, istate);
						return;
					}
					continue;
				}
				/*
				   itask = 4.
				   See if tout or tcrit was reached.  Adjust h if necessary.
				 */
				if (itask == 4) {
					if ((tn - tout) * h >= 0.) {
						intdy(tout, 0, y, &iflag);
						*t = tout;
						*istate = 2;
						illin = 0;
						return;
					} else {
						hmx = fabs(tn) + fabs(h);
						ihit = fabs(tn - tcrit) <= (100. * ETA * hmx);
						if (ihit) {
							successreturn(y, t, itask, ihit, tcrit, istate);
							return;
						}
						tnext = tn + h * (1. + 4. * ETA);
						if ((tnext - tcrit) * h <= 0.)
							continue;
						h = (tcrit - tn) * (1. - 4. * ETA);
						jstart = -2;
						continue;
					}
				}	/* end if ( itask == 4 )   */
				/*
				   itask = 5.
				   See if tcrit was reached and jump to exit.
				 */
				if (itask == 5) {
					hmx = fabs(tn) + fabs(h);
					ihit = fabs(tn - tcrit) <= (100. * ETA * hmx);
					successreturn(y, t, itask, ihit, tcrit, istate);
					return;
				}
			}		/* end if ( kflag == 0 )   */
			/*
			   kflag = -1, error test failed repeatedly or with fabs(h) = hmin.
			   kflag = -2, convergence failed repeatedly or with fabs(h) = hmin.
			 */
			if (kflag == -1 || kflag == -2) {
				fprintf(stderr, "lsoda -- at t = %g and step size h = %g, the\n", tn, h);
				if (kflag == -1) {
					fprintf(stderr, "         error test failed repeatedly or\n");
					fprintf(stderr, "         with fabs(h) = hmin\n");
					*istate = -4;
				}
				if (kflag == -2) {
					fprintf(stderr, "         corrector convergence failed repeatedly or\n");
					fprintf(stderr, "         with fabs(h) = hmin\n");
					*istate = -5;
				}
				big = 0.;
				imxer = 1;
				for (i = 1; i <= n; i++) {
					size = fabs(vec.acor[i]) * vec.ewt[i];
					if (big < size) {
						big = size;
						imxer = i;
					}
				}
				terminate2(y, t);
				return;
			}		/* end if ( kflag == -1 || kflag == -2 )   */
		}			/* end while   */

	}				/* end lsoda   */

