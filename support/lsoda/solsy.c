#include "common.h"
#include <stdio.h>
#include "blas.h"

int solsy(struct common_t * common, int neq, double *y)

/*
   This routine manages the solution of the linear system arising from
   a chord iteration.  It is called if _C(miter) != 0.
   If _C(miter) is 2, it calls dgesl to accomplish this.
   If _C(miter) is 5, it calls dgbsl.

   y = the right-hand side vector on input, and the solution vector
       on output.
*/

{
	if (_C(miter) != 2) {
		printf("solsy -- _C(miter) != 2\n");
		return 0;
	}
	if (_C(miter) == 2)
		dgesl(_C(wm), neq, _C(ipvt), y, 0);
	return 1;

}

