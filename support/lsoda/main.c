#include <stdio.h>
#include "lsoda.h"

static void fex(int * neq, double *t, double *y, double *ydot)
{
	ydot[0] = 1.0E4 * y[1] * y[2] - .04E0 * y[0];
	ydot[2] = 3.0E7 * y[1] * y[1];
	ydot[1] = -1.0 * (ydot[0] + ydot[2]);
//	printf("fex called, t = %g neq = %d y = %g %g %g ydot= %g %g %g \n", *t , *neq,  y[0], y[1], y[2], ydot[0], ydot[1], ydot[2]);
}

int main() {
	double          atol[3], rtol[3], t, tout, y[3], yout[3];
	int             neq = 3;
	int             iout, istate;

	y[0] = 1.0E0;
	y[1] = 0.0E0;
	y[2] = 0.0E0;
	
	t = 0.0E0;
	tout = 0.4E0;
	rtol[0] = rtol[2] = 1.0E-4;
	rtol[1] = 1.0E-8;
	atol[0] = 1.0E-6;
	atol[1] = 1.0E-10;
	atol[2] = 1.0E-6;

	for (iout = 1; iout <= 12; iout++) {
		istate = lsoda(y, yout, neq, &t, tout, fex, rtol, 1, atol, 3, NULL, NULL, NULL, NULL);
		printf(" at t= %12.4e y= %14.6e %14.6e %14.6e\n", t, yout[0], yout[1], yout[2]);
		if (istate <= 0) {
			printf("error istate = %d\n", istate);
			exit(0);
		}
		y[0] = yout[0];
		y[1] = yout[1];
		y[2] = yout[2];
		tout = tout * 10.0E0;
	}

	return 0;
}
