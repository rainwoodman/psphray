#include "lsoda.h"
#include "lsoda_internal.h"
#include "common.h"
#include <math.h>
void scaleh(int neq, double *rh, double *pdh, double hmxi)
{
	double          r;
	int             j, i;
	double ** yh = vec.yh;
/*
   If h is being changed, the h ratio rh is checked against rmax, hmin,
   and hmxi, and the yh array is rescaled.  ialth is set to (nq + 1) = nq + 1
   to prevent a change of h for that many steps, unless forced by a
   convergence or error test failure.
*/
	*rh = min(*rh, rmax);
	*rh = *rh / max(1., fabs(h) * hmxi * *rh);
/*
   If meth = 1, also restrict the new step size by the stability region.
   If this reduces h, set irflag to 1 so that if there are roundoff
   problems later, we can assume that is the cause of the trouble.
*/
	if (meth == 1) {
		irflag = 0;
		*pdh = max(fabs(h) * pdlast, 0.000001);
		if ((*rh * *pdh * 1.00001) >= sm1[nq]) {
			*rh = sm1[nq] / *pdh;
			irflag = 1;
		}
	}
	r = 1.;
	for (j = 2; j <= (nq + 1); j++) {
		r *= *rh;
		for (i = 1; i <= neq; i++)
			yh[j][i] *= r;
	}
	h *= *rh;
	rc *= *rh;
	ialth = (nq + 1);

}				/* end scaleh   */
