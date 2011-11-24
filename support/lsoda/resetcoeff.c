#include "common.h"

void resetcoeff()
/*
   The el vector and related constants are reset
   whenever the order nq is changed, or at the start of the problem.
*/
{
	int             i;
	double         *ep1;

	ep1 = elco[nq];
	for (i = 1; i <= l; i++)
		el[i] = ep1[i];
	rc = rc * el[1] / el0;
	el0 = el[1];
	conit = 0.5 / (double) (nq + 2);

}

