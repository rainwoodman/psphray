
#include "common.h"
#include <math.h>
void ewset(int itol, double *rtol, double *atol, double *ycur)
{
	int             i;

	switch (itol) {
	case 1:
		for (i = 1; i <= n; i++)
			ewt[i] = rtol[1] * fabs(ycur[i]) + atol[1];
		break;
	case 2:
		for (i = 1; i <= n; i++)
			ewt[i] = rtol[1] * fabs(ycur[i]) + atol[i];
		break;
	case 3:
		for (i = 1; i <= n; i++)
			ewt[i] = rtol[i] * fabs(ycur[i]) + atol[1];
		break;
	case 4:
		for (i = 1; i <= n; i++)
			ewt[i] = rtol[i] * fabs(ycur[i]) + atol[i];
		break;
	}

}				/* end ewset   */

