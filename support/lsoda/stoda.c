#include "lsoda.h"
#include "common.h"
#include "lsoda_internal.h"

#include <math.h>
#include "blas.h"

/*
   This routine returns from stoda to lsoda.  Hence freevectors() is
   not executed.
*/
void endstoda()
{
	double          r;
	int             i;

	r = 1. / tesco[nqu][2];
	for (i = 1; i <= n; i++)
		acor[i] *= r;
	hold = h;

}

int stoda(int neq, double *y, _lsoda_f f, void *_data, int jstart, double hmxi, double hmin, int mxords, int mxordn)
{
	int kflag;
	int             corflag, orderflag;
	int             i, i1, j, m, ncf;
	double          del, delp, dsm, dup, exup, r, rh, rhup, told;
	double          pdh, pnorm;

/*
   stoda performs one step of the integration of an initial value
   problem for a system of ordinary differential equations.
   Note.. stoda is independent of the value of the iteration method
   indicator miter, when this is != 0, and hence is independent
   of the type of chord method used, or the Jacobian structure.
   Communication with stoda is done with the following variables:

   jstart = an integer used for input only, with the following
            values and meanings:

               0  perform the first step,
             > 0  take a new step continuing from the last,
              -1  take the next step with a new value of h,
                  n, meth, miter, and/or matrix parameters.
              -2  take the next step with a new value of h,
                  but with other inputs unchanged.

   kflag = a completion code with the following meanings:

             0  the step was successful,
            -1  the requested error could not be achieved,
            -2  corrector convergence could not be achieved,
            -3  fatal error in prja or solsy.

   miter = corrector iteration method:

             0  functional iteration,
            >0  a chord method corresponding to jacobian type jt.

*/
	kflag = 0;
	told = tn;
	ncf = 0;
	ierpj = 0;
	iersl = 0;
	jcur = 0;
	delp = 0.;

/*
   On the first call, the order is set to 1, and other variables are
   initialized.  rmax is the maximum ratio by which h can be increased
   in a single step.  It is initially 1.e4 to compensate for the small
   initial h, but then is normally equal to 10.  If a filure occurs
   (in corrector convergence or error test), rmax is set at 2 for
   the next increase.
   cfode is called to get the needed coefficients for both methods.
*/
	if (jstart == 0) {
		lmax = maxord + 1;
		nq = 1;
		l = 2;
		ialth = 2;
		rmax = 10000.;
		rc = 0.;
		el0 = 1.;
		crate = 0.7;
		hold = h;
		nslp = 0;
		ipup = miter;
/*
   Initialize switching parameters.  meth = 1 is assumed initially.
*/
		icount = 20;
		irflag = 0;
		pdest = 0.;
		pdlast = 0.;
		ratio = 5.;
		cfode(2);
		for (i = 1; i <= 5; i++)
			cm2[i] = tesco[i][2] * elco[i][i + 1];
		cfode(1);
		for (i = 1; i <= 12; i++)
			cm1[i] = tesco[i][2] * elco[i][i + 1];
		resetcoeff();
	}			/* end if ( jstart == 0 )   */
	/*
	   The following block handles preliminaries needed when jstart = -1.
	   ipup is set to miter to force a matrix update.
	   If an order increase is about to be considered ( ialth = 1 ),
	   ialth is reset to 2 to postpone consideration one more step.
	   If the caller has changed meth, cfode is called to reset
	   the coefficients of the method.
	   If h is to be changed, yh must be rescaled.
	   If h or meth is being changed, ialth is reset to l = nq + 1
	   to prevent further changes in h for that many steps.
	*/
	if (jstart == -1) {
		ipup = miter;
		lmax = maxord + 1;
		if (ialth == 1)
			ialth = 2;
		if (meth != mused) {
			cfode(meth);
			ialth = l;
			resetcoeff();
		}
		if (h != hold) {
			rh = h / hold;
			h = hold;
			scaleh(&rh, &pdh, hmxi);
		}
	}			/* if ( jstart == -1 )   */
	if (jstart == -2) {
		if (h != hold) {
			rh = h / hold;
			h = hold;
			scaleh(&rh, &pdh, hmxi);
		}
	}			/* if ( jstart == -2 )   */
	/*
	   Prediction.
	   This section computes the predicted values by effectively
	   multiplying the yh array by the pascal triangle matrix.
	   rc is the ratio of new to old values of the coefficient h * el[1].
	   When rc differs from 1 by more than ccmax, ipup is set to miter
	   to force pjac to be called, if a jacobian is involved.
	   In any case, prja is called at least every msbp steps.
	*/
	while (1) {
		while (1) {
			if (fabs(rc - 1.) > ccmax)
				ipup = miter;
			if (nst >= nslp + msbp)
				ipup = miter;
			tn += h;
			for (j = nq; j >= 1; j--)
				for (i1 = j; i1 <= nq; i1++) {
					yp1 = yh[i1];
					yp2 = yh[i1 + 1];
					for (i = 1; i <= n; i++)
						yp1[i] += yp2[i];
				}
			pnorm = vmnorm(n, yh[1], ewt);

			correction(neq, y, f, &corflag, pnorm, &del, &delp, &told, &ncf, &rh, &m, hmin, _data);
			if (corflag == 0)
				break;
			if (corflag == 1) {
				rh = max(rh, hmin / fabs(h));
				scaleh(&rh, &pdh, hmxi);
				continue;
			}
			if (corflag == 2) {
				kflag = -2;
				hold = h;
				jstart = 1;
				return kflag;
			}
		}		/* end inner while ( corrector loop )   */
/*
   The corrector has converged.  jcur is set to 0
   to signal that the Jacobian involved may need updating later.
   The local error test is done now.
*/
		jcur = 0;
		if (m == 0)
			dsm = del / tesco[nq][2];
		if (m > 0)
			dsm = vmnorm(n, acor, ewt) / tesco[nq][2];
		if (dsm <= 1.) {
/*
   After a successful step, update the yh array.
   Decrease icount by 1, and if it is -1, consider switching methods.
   If a method switch is made, reset various parameters,
   rescale the yh array, and exit.  If there is no switch,
   consider changing h if ialth = 1.  Otherwise decrease ialth by 1.
   If ialth is then 1 and nq < maxord, then acor is saved for
   use in a possible order increase on the next step.
   If a change in h is considered, an increase or decrease in order
   by one is considered also.  A change in h is made only if it is by
   a factor of at least 1.1.  If not, ialth is set to 3 to prevent
   testing for that many steps.
*/
			kflag = 0;
			nst++;
			hu = h;
			nqu = nq;
			mused = meth;
			for (j = 1; j <= l; j++) {
				yp1 = yh[j];
				r = el[j];
				for (i = 1; i <= n; i++)
					yp1[i] += r * acor[i];
			}
			icount--;
			if (icount < 0) {
				methodswitch(dsm, pnorm, &pdh, &rh, mxords, mxordn);
				if (meth != mused) {
					rh = max(rh, hmin / fabs(h));
					scaleh(&rh, &pdh, hmxi);
					rmax = 10.;
					endstoda();
					break;
				}
			}
/*
   No method switch is being made.  Do the usual step/order selection.
*/
			ialth--;
			if (ialth == 0) {
				rhup = 0.;
				if (l != lmax) {
					yp1 = yh[lmax];
					for (i = 1; i <= n; i++)
						savf[i] = acor[i] - yp1[i];
					dup = vmnorm(n, savf, ewt) / tesco[nq][3];
					exup = 1. / (double) (l + 1);
					rhup = 1. / (1.4 * pow(dup, exup) + 0.0000014);
				}
				orderswitch(&rhup, dsm, &pdh, &rh, &orderflag, kflag);
/*
   No change in h or nq.
*/
				if (orderflag == 0) {
					endstoda();
					break;
				}
/*
   h is changed, but not nq.
*/
				if (orderflag == 1) {
					rh = max(rh, hmin / fabs(h));
					scaleh(&rh, &pdh, hmxi);
					rmax = 10.;
					endstoda();
					break;
				}
/*
   both nq and h are changed.
*/
				if (orderflag == 2) {
					resetcoeff();
					rh = max(rh, hmin / fabs(h));
					scaleh(&rh, &pdh, hmxi);
					rmax = 10.;
					endstoda();
					break;
				}
			}	/* end if ( ialth == 0 )   */
			if (ialth > 1 || l == lmax) {
				endstoda();
				break;
			}
			yp1 = yh[lmax];
			for (i = 1; i <= n; i++)
				yp1[i] = acor[i];
			endstoda();
			break;
		}
		/* end if ( dsm <= 1. )   */
		/*
		   The error test failed.  kflag keeps track of multiple failures.
		   Restore tn and the yh array to their previous values, and prepare
		   to try the step again.  Compute the optimum step size for this or
		   one lower.  After 2 or more failures, h is forced to decrease
		   by a factor of 0.2 or less.
		 */ 
		else {
			kflag--;
			tn = told;
			for (j = nq; j >= 1; j--)
				for (i1 = j; i1 <= nq; i1++) {
					yp1 = yh[i1];
					yp2 = yh[i1 + 1];
					for (i = 1; i <= n; i++)
						yp1[i] -= yp2[i];
				}
			rmax = 2.;
			if (fabs(h) <= hmin * 1.00001) {
				kflag = -1;
				hold = h;
				jstart = 1;
				break;
			}
			if (kflag > -3) {
				rhup = 0.;
				orderswitch(&rhup, dsm, &pdh, &rh, &orderflag, kflag);
				if (orderflag == 1 || orderflag == 0) {
					if (orderflag == 0)
						rh = min(rh, 0.2);
					rh = max(rh, hmin / fabs(h));
					scaleh(&rh, &pdh, hmxi);
				}
				if (orderflag == 2) {
					resetcoeff();
					rh = max(rh, hmin / fabs(h));
					scaleh(&rh, &pdh, hmxi);
				}
				continue;
			}
			/* if ( kflag > -3 )   */
			/*
			   Control reaches this section if 3 or more failures have occurred.
			   If 10 failures have occurred, exit with kflag = -1.
			   It is assumed that the derivatives that have accumulated in the
			   yh array have errors of the wrong order.  Hence the first
			   derivative is recomputed, and the order is set to 1.  Then
			   h is reduced by a factor of 10, and the step is retried,
			   until it succeeds or h reaches hmin.
			 */ 
			else {
				if (kflag == -10) {
					kflag = -1;
					hold = h;
					jstart = 1;
					break;
				} else {
					rh = 0.1;
					rh = max(hmin / fabs(h), rh);
					h *= rh;
					yp1 = yh[1];
					for (i = 1; i <= n; i++)
						y[i] = yp1[i];
					(*f) (tn, y + 1, savf + 1, _data);
					nfe++;
					yp1 = yh[2];
					for (i = 1; i <= n; i++)
						yp1[i] = h * savf[i];
					ipup = miter;
					ialth = 5;
					if (nq == 1)
						continue;
					nq = 1;
					l = 2;
					resetcoeff();
					continue;
				}
			}	/* end else -- kflag <= -3 */
		}		/* end error failure handling   */
	}			/* end outer while   */

	return kflag;
}				/* end stoda   */


