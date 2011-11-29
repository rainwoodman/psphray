#include "common.h"

void resetcoeff(struct common_t * common)
/*
   The _C(el) vector and related constants are reset
   whenever the order _C(nq) is changed, or at the start of the problem.
*/
{
	int             i;

	double el0 = _C(el)[1];
	for (i = 1; i <= (_C(nq) + 1); i++)
		_C(el)[i] = _C(elco)[_C(nq)][i];
	_C(rc) = _C(rc) * _C(el)[1] / el0;
}

