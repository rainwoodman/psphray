#include "common.h"
#include <stdio.h>
#include "blas.h"

void solsy(int neq, double *y, double ** wm, int * ipvt)

/*
   This routine manages the solution of the linear system arising from
   a chord iteration.  It is called if miter != 0.
   If miter is 2, it calls dgesl to accomplish this.
   If miter is 5, it calls dgbsl.

   y = the right-hand side vector on input, and the solution vector
       on output.
*/

{
	iersl = 0;
	if (miter != 2) {
		printf("solsy -- miter != 2\n");
		return;
	}
	if (miter == 2)
		dgesl(wm, neq, ipvt, y, 0);
	return;

}

