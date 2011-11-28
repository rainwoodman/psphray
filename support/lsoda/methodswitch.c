#include <math.h>
#include "lsoda.h"
#include "lsoda_internal.h"
#include "common.h"
#include "blas.h"

void methodswitch(int neq, double dsm, double pnorm, double *pdh, double *rh, int mxords, int mxordn)
{
	int             lm1, lm1p1, lm2, lm2p1, nqm1, nqm2;
	double          rh1, rh2, rh1it, exm2, dm2, exm1, dm1, alpha, exsm;
	double ** yh = vec.yh;
	double * ewt = vec.ewt;
/*
   We are current using an Adams method.  Consider switching to bdf.
   If the current order is greater than 5, assume the problem is
   not stiff, and skip this section.
   If the Lipschitz constant and error estimate are not polluted
   by roundoff, perform the usual test.
   Otherwise, switch to the bdf methods if the last step was
   restricted to insure stability ( irflag = 1 ), and stay with Adams
   method if not.  When switching to bdf with polluted error estimates,
   in the absence of other information, double the step size.

   When the estimates are ok, we make the usual test by computing
   the step size we could have (ideally) used on this step,
   with the current (Adams) method, and also that for the bdf.
   If nq > mxords, we consider changing to order mxords on switching.
   Compare the two step sizes to decide whether to switch.
   The step size advantage must be at least ratio = 5 to switch.
*/
	if (meth == 1) {
		if (nq > 5)
			return;
		if (dsm <= (100. * pnorm * ETA) || pdest == 0.) {
			if (irflag == 0)
				return;
			rh2 = 2.;
			nqm2 = min(nq, mxords);
		} else {
			exsm = 1. / (double) (nq + 1);
			rh1 = 1. / (1.2 * pow(dsm, exsm) + 0.0000012);
			rh1it = 2. * rh1;
			*pdh = pdlast * fabs(h);
			if ((*pdh * rh1) > 0.00001)
				rh1it = sm1[nq] / *pdh;
			rh1 = min(rh1, rh1it);
			if (nq > mxords) {
				nqm2 = mxords;
				lm2 = mxords + 1;
				exm2 = 1. / (double) lm2;
				lm2p1 = lm2 + 1;
				dm2 = vmnorm(neq, yh[lm2p1], ewt) / cm2[mxords];
				rh2 = 1. / (1.2 * pow(dm2, exm2) + 0.0000012);
			} else {
				dm2 = dsm * (cm1[nq] / cm2[nq]);
				rh2 = 1. / (1.2 * pow(dm2, exsm) + 0.0000012);
				nqm2 = nq;
			}
			if (rh2 < ratio * rh1)
				return;
		}
/*
   The method switch test passed.  Reset relevant quantities for bdf.
*/
		*rh = rh2;
		icount = 20;
		meth = 2;
		miter = 2;
		pdlast = 0.;
		nq = nqm2;
		return;
	}			/* end if ( meth == 1 )   */
	/*
	   We are currently using a bdf method, considering switching to Adams.
	   Compute the step size we could have (ideally) used on this step,
	   with the current (bdf) method, and also that for the Adams.
	   If nq > mxordn, we consider changing to order mxordn on switching.
	   Compare the two step sizes to decide whether to switch.
	   The step size advantage must be at least 5/ratio = 1 to switch.
	   If the step size for Adams would be so small as to cause
	   roundoff pollution, we stay with bdf.
	*/
	exsm = 1. / (double) (nq + 1);
	if (mxordn < nq) {
		nqm1 = mxordn;
		lm1 = mxordn + 1;
		exm1 = 1. / (double) lm1;
		lm1p1 = lm1 + 1;
		dm1 = vmnorm(neq, yh[lm1p1], ewt) / cm1[mxordn];
		rh1 = 1. / (1.2 * pow(dm1, exm1) + 0.0000012);
	} else {
		dm1 = dsm * (cm2[nq] / cm1[nq]);
		rh1 = 1. / (1.2 * pow(dm1, exsm) + 0.0000012);
		nqm1 = nq;
		exm1 = exsm;
	}
	rh1it = 2. * rh1;
	*pdh = pdnorm * fabs(h);
	if ((*pdh * rh1) > 0.00001)
		rh1it = sm1[nqm1] / *pdh;
	rh1 = min(rh1, rh1it);
	rh2 = 1. / (1.2 * pow(dsm, exsm) + 0.0000012);
	if ((rh1 * ratio) < (5. * rh2))
		return;
	alpha = max(0.001, rh1);
	dm1 *= pow(alpha, exm1);
	if (dm1 <= 1000. * ETA * pnorm)
		return;
/*
   The switch test passed.  Reset relevant quantities for Adams.
*/
	*rh = rh1;
	icount = 20;
	meth = 1;
	miter = 0;
	pdlast = 0.;
	nq = nqm1;

}				/* end methodswitch   */

