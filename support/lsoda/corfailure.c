#include "lsoda.h"
#include "common.h"
#include "lsoda_internal.h"
#include <math.h>

int corfailure(int neq, double *told, double *rh, int *ncf, double hmin)
{
	int             j, i1, i;
	double ** yh = vec.yh;

	ncf++;
	rmax = 2.;
	tn = *told;
	for (j = nq; j >= 1; j--)
		for (i1 = j; i1 <= nq; i1++) {
			for (i = 1; i <= neq; i++)
				yh[i1][i] -= yh[i1 + 1][i];
		}
	if (fabs(h) <= hmin * 1.00001 || *ncf == mxncf) {
		return 2;
	}
	*rh = 0.25;
	ipup = miter;
	return 1;
}

