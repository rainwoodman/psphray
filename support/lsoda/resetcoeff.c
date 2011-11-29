#include "common.h"

void resetcoeff()
/*
   The el vector and related constants are reset
   whenever the order nq is changed, or at the start of the problem.
*/
{
	int             i;

	double el0 = el[1];
	for (i = 1; i <= (nq + 1); i++)
		el[i] = elco[nq][i];
	rc = rc * el[1] / el0;
}

