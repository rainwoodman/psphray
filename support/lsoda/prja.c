#include <stdio.h>
#include <math.h>
#include "lsoda.h"
#include "lsoda_internal.h"
#include "blas.h"
#include "common.h"

void prja(int neq, double *y, _lsoda_f f, void *_data)
{
	int             i, ier, j;
	double          fac, hl0, r, r0, yj;
	double * savf = vec.savf;
	double * ewt = vec.ewt;
	double * acor = vec.acor;
	double ** wm = vec.wm;
	int * ipvt = vec.ipvt;

/*
   prja is called by stoda to compute and process the matrix
   P = I - h * el[1] * J, where J is an approximation to the Jacobian.
   Here J is computed by finite differencing.
   J, scaled by -h * el[1], is stored in wm.  Then the norm of J ( the
   matrix norm consistent with the weighted max-norm on vectors given
   by vmnorm ) is computed, and J is overwritten by P.  P is then
   subjected to LU decomposition in preparation for later solution
   of linear systems with p as coefficient matrix.  This is done
   by dgefa if miter = 2, and by dgbfa if miter = 5.
*/
	nje++;
	ierpj = 0;
	jcur = 1;
	hl0 = h * el0;
/*
   If miter = 2, make n calls to f to approximate J.
*/
	if (miter != 2) {
		fprintf(stderr, "[prja] miter != 2\n");
		return;
	}
	if (miter == 2) {
		fac = vmnorm(n, savf, ewt);
		r0 = 1000. * fabs(h) * ETA * ((double) n) * fac;
		if (r0 == 0.)
			r0 = 1.;
		for (j = 1; j <= n; j++) {
			yj = y[j];
			r = fmax(SQRTETA * fabs(yj), r0 / ewt[j]);
			y[j] += r;
			fac = -hl0 / r;
			(*f) (tn, y + 1, acor + 1, _data);
			for (i = 1; i <= n; i++)
				wm[i][j] = (acor[i] - savf[i]) * fac;
			y[j] = yj;
		}
		nfe += n;
/*
   Compute norm of Jacobian.
*/
		pdnorm = fnorm(n, wm, ewt) / fabs(hl0);
/*
   Add identity matrix.
*/
		for (i = 1; i <= n; i++)
			wm[i][i] += 1.;
/*
   Do LU decomposition on P.
*/
		dgefa(wm, n, ipvt, &ier);
		if (ier != 0)
			ierpj = 1;
		return;
	}
}				/* end prja   */

