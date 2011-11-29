#include <stdio.h>
#include <math.h>
#include "lsoda.h"
#include "lsoda_internal.h"
#include "common.h"

void intdy(int neq, double t, int k, double *dky, int *iflag)

/*
   Intdy computes interpolated values of the k-th derivative of the
   dependent variable vector y, and stores it in dky.  This routine
   is called within the package with k = 0 and *t = tout, but may
   also be called by the user for any k up to the current order.
   ( See detailed instructions in the usage documentation. )

   The computed values in dky are gotten by interpolation using the
   Nordsieck history array yh.  This array corresponds uniquely to a
   vector-valued polynomial of degree nqcur or less, and dky is set
   to the k-th derivative of this polynomial at t.
   The formula for dky is

             q
   dky[i] = sum c[k][j] * ( t - tn )^(j-k) * h^(-j) * yh[j+1][i]
            j=k

   where c[k][j] = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
   The quantities nq = nqcur, l = nq+1, neq = neq, tn, and h are declared
   static globally.  The above sum is done in reverse order.
   *iflag is returned negative if either k or t is out of bounds.
*/

{
	int             i, ic, j, jj, jp1;
	double          c, r, s, tp;

	*iflag = 0;
	if (k < 0 || k > nq) {
		fprintf(stderr, "[intdy] k = %d illegal\neq", k);
		*iflag = -1;
		return;
	}
	tp = tn - hu - 100. * ETA * (tn + hu);
	if ((t - tp) * (t - tn) > 0.) {
		fprintf(stderr, "intdy -- t = %g illegal. t not in interval tcur - hu to tcur\neq", t);
		*iflag = -2;
		return;
	}
	s = (t - tn) / h;
	ic = 1;
	for (jj = (nq + 1) - k; jj <= nq; jj++)
		ic *= jj;
	c = (double) ic;
	for (i = 1; i <= neq; i++)
		dky[i] = c * yh[nq + 1][i];
	for (j = nq - 1; j >= k; j--) {
		jp1 = j + 1;
		ic = 1;
		for (jj = jp1 - k; jj <= j; jj++)
			ic *= jj;
		c = (double) ic;
		for (i = 1; i <= neq; i++)
			dky[i] = c * yh[jp1][i] + s * dky[i];
	}
	if (k == 0)
		return;
	r = pow(h, (double) (-k));
	for (i = 1; i <= neq; i++)
		dky[i] *= r;

}				/* end intdy   */

