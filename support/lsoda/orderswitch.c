#include "lsoda.h"
#include "lsoda_internal.h"
#include "common.h"
#include <math.h>
#include "blas.h"

void orderswitch(double *rhup, double dsm, double *pdh, double *rh, int *orderflag, int kflag)

/*
   Regardless of the success or failure of the step, factors
   rhdn, rhsm, and rhup are computed, by which h could be multiplied
   at order nq - 1, order nq, or order nq + 1, respectively.
   In the case of a failure, rhup = 0. to avoid an order increase.
   The largest of these is determined and the new order chosen
   accordingly.  If the order is to be increased, we compute one
   additional scaled derivative.

   orderflag = 0  : no change in h or nq,
               1  : change in h but not nq,
               2  : change in both h and nq.
*/

{
	int             newq, i;
	double          exsm, rhdn, rhsm, ddn, exdn, r;
	double ** yh = vec.yh;
	double * ewt = vec.ewt;
	double * acor = vec.acor;

	*orderflag = 0;

	exsm = 1. / (double) l;
	rhsm = 1. / (1.2 * pow(dsm, exsm) + 0.0000012);

	rhdn = 0.;
	if (nq != 1) {
		ddn = vmnorm(n, yh[l], ewt) / tesco[nq][1];
		exdn = 1. / (double) nq;
		rhdn = 1. / (1.3 * pow(ddn, exdn) + 0.0000013);
	}
/*
   If meth = 1, limit rh accordinfg to the stability region also.
*/
	if (meth == 1) {
		*pdh = max(fabs(h) * pdlast, 0.000001);
		if (l < lmax)
			*rhup = min(*rhup, sm1[l] / *pdh);
		rhsm = min(rhsm, sm1[nq] / *pdh);
		if (nq > 1)
			rhdn = min(rhdn, sm1[nq - 1] / *pdh);
		pdest = 0.;
	}
	if (rhsm >= *rhup) {
		if (rhsm >= rhdn) {
			newq = nq;
			*rh = rhsm;
		} else {
			newq = nq - 1;
			*rh = rhdn;
			if (kflag < 0 && *rh > 1.)
				*rh = 1.;
		}
	} else {
		if (*rhup <= rhdn) {
			newq = nq - 1;
			*rh = rhdn;
			if (kflag < 0 && *rh > 1.)
				*rh = 1.;
		} else {
			*rh = *rhup;
			if (*rh >= 1.1) {
				r = el[l] / (double) l;
				nq = l;
				l = nq + 1;
				yp1 = yh[l];
				for (i = 1; i <= n; i++)
					yp1[i] = acor[i] * r;
				*orderflag = 2;
				return;
			} else {
				ialth = 3;
				return;
			}
		}
	}
/*
   If meth = 1 and h is restricted by stability, bypass 10 percent test.
*/
	if (meth == 1) {
		if ((*rh * *pdh * 1.00001) < sm1[newq])
			if (kflag == 0 && *rh < 1.1) {
				ialth = 3;
				return;
			}
	} else {
		if (kflag == 0 && *rh < 1.1) {
			ialth = 3;
			return;
		}
	}
	if (kflag <= -2)
		*rh = min(*rh, 0.2);
/*
   If there is a change of order, reset nq, l, and the coefficients.
   In any case h is reset according to rh and the yh array is rescaled.
   Then exit or redo the step.
*/
	if (newq == nq) {
		*orderflag = 1;
		return;
	}
	nq = newq;
	l = nq + 1;
	*orderflag = 2;

}				/* end orderswitch   */


