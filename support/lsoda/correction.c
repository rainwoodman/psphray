
#include <math.h>
#include "common.h"
#include "blas.h"
#include "lsoda_internal.h"

void correction(int neq, double *y, _lsoda_f f, int *corflag, double pnorm, double *del, double *delp, double *told,
					   int *ncf, double *rh, int *m, void *_data)
/*
   *corflag = 0 : corrector converged,
              1 : step size to be reduced, redo prediction,
              2 : corrector cannot converge, failure flag.
*/

{
	int             i;
	double          rm, rate, dcon;

/*
   Up to maxcor corrector iterations are taken.  A convergence test is
   made on the r.m.s. norm of each correction, weighted by the error
   weight vector ewt.  The sum of the corrections is accumulated in the
   vector acor[i].  The yh array is not altered in the corrector loop.
*/

	*m = 0;
	*corflag = 0;
	rate = 0.;
	*del = 0.;
	yp1 = yh[1];
	for (i = 1; i <= n; i++)
		y[i] = yp1[i];
	(*f) (tn, y + 1, savf + 1, _data);
	nfe++;
/*
   If indicated, the matrix P = I - h * el[1] * J is reevaluated and
   preprocessed before starting the corrector iteration.  ipup is set
   to 0 as an indicator that this has been done.
*/
	while (1) {
		if (*m == 0) {
			if (ipup > 0) {
				prja(neq, y, f, _data);
				ipup = 0;
				rc = 1.;
				nslp = nst;
				crate = 0.7;
				if (ierpj != 0) {
					corfailure(told, rh, ncf, corflag);
					return;
				}
			}
			for (i = 1; i <= n; i++)
				acor[i] = 0.;
		}		/* end if ( *m == 0 )   */
		if (miter == 0) {
/*
   In case of functional iteration, update y directly from
   the result of the last function evaluation.
*/
			yp1 = yh[2];
			for (i = 1; i <= n; i++) {
				savf[i] = h * savf[i] - yp1[i];
				y[i] = savf[i] - acor[i];
			}
			*del = vmnorm(n, y, ewt);
			yp1 = yh[1];
			for (i = 1; i <= n; i++) {
				y[i] = yp1[i] + el[1] * savf[i];
				acor[i] = savf[i];
			}
		}
		/* end functional iteration   */
		/*
		   In the case of the chord method, compute the corrector error,
		   and solve the linear system with that as right-hand side and
		   P as coefficient matrix.
		 */ 
		else {
			yp1 = yh[2];
			for (i = 1; i <= n; i++)
				y[i] = h * savf[i] - (yp1[i] + acor[i]);
			solsy(y);
			*del = vmnorm(n, y, ewt);
			yp1 = yh[1];
			for (i = 1; i <= n; i++) {
				acor[i] += y[i];
				y[i] = yp1[i] + el[1] * acor[i];
			}
		}		/* end chord method   */
/*
   Test for convergence.  If *m > 0, an estimate of the convergence
   rate constant is stored in crate, and this is used in the test.

   We first check for a change of iterates that is the size of
   roundoff error.  If this occurs, the iteration has converged, and a
   new rate estimate is not formed.
   In all other cases, force at least two iterations to estimate a
   local Lipschitz constant estimate for Adams method.
   On convergence, form pdest = local maximum Lipschitz constant
   estimate.  pdlast is the most recent nonzero estimate.
*/
		if (*del <= 100. * pnorm * ETA)
			break;
		if (*m != 0 || meth != 1) {
			if (*m != 0) {
				rm = 1024.0;
				if (*del <= (1024. * *delp))
					rm = *del / *delp;
				rate = max(rate, rm);
				crate = max(0.2 * crate, rm);
			}
			dcon = *del * min(1., 1.5 * crate) / (tesco[nq][2] * conit);
			if (dcon <= 1.) {
				pdest = max(pdest, rate / fabs(h * el[1]));
				if (pdest != 0.)
					pdlast = pdest;
				break;
			}
		}
/*
   The corrector iteration failed to converge.
   If miter != 0 and the Jacobian is out of date, prja is called for
   the next try.   Otherwise the yh array is retracted to its values
   before prediction, and h is reduced, if possible.  If h cannot be
   reduced or mxncf failures have occured, exit with corflag = 2.
*/
		(*m)++;
		if (*m == maxcor || (*m >= 2 && *del > 2. * *delp)) {
			if (miter == 0 || jcur == 1) {
				corfailure(told, rh, ncf, corflag);
				return;
			}
			ipup = miter;
/*
   Restart corrector if Jacobian is recomputed.
*/
			*m = 0;
			rate = 0.;
			*del = 0.;
			yp1 = yh[1];
			for (i = 1; i <= n; i++)
				y[i] = yp1[i];
			(*f) (tn, y + 1, savf + 1, _data);
			nfe++;
		}
/*
   Iterate corrector.
*/
		else {
			*delp = *del;
			(*f) (tn, y + 1, savf + 1, _data);
			nfe++;
		}
	}			/* end while   */
}				/* end correction   */

