
#include <math.h>
#include <stdio.h>
int ewset(struct common_t * common, const int neq, double * ewt, const double *rtol, const double *atol, const double *ycur)
{
	int             i;

	for (i = 1; i <= neq; i++)
		ewt[i] = rtol[i] * fabs(ycur[i]) + atol[i];

	for (i = 1; i <= neq; i++) {
/*
		if (ewt[i] <= 0.) {
			fprintf(stderr, "[lsoda] ewt[%d] = %g <= 0.\neq", i, ewt[i]);
			return 0;
		}
*/
		ewt[i] = 1. / ewt[i];
	}
	return 1;
}				/* end ewset   */

