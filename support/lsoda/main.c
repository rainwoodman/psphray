#include <stdio.h>
#include "common.h"
#include "lsoda.h"

static void fex(double t, double *y, double *ydot, void *data)
{
	ydot[0] = 1.0E4 * y[1] * y[2] - .04E0 * y[0];
	ydot[2] = 3.0E7 * y[1] * y[1];
	ydot[1] = -1.0 * (ydot[0] + ydot[2]);
}

int main(void)
{
	double          atol[4], rtol[4], t, tout, y[4];
	int             neq = 3;
	int             itol, itask, istate, iopt, jt, iout;

	y[1] = 1.0E0;
	y[2] = 0.0E0;
	y[3] = 0.0E0;
	t = 0.0E0;
	tout = 0.4E0;
	itol = 2;
	rtol[0] = 0.0;
	atol[0] = 0.0;
	rtol[1] = rtol[3] = 1.0E-4;
	rtol[2] = 1.0E-8;
	atol[1] = 1.0E-6;
	atol[2] = 1.0E-10;
	atol[3] = 1.0E-6;
	itask = 1;
	istate = 1;
	iopt = 0;
	jt = 2;

	struct lsoda_opt_t opt = {0};
	struct lsoda_context_t ctx = {
		.function = fex,
		.neq = neq,
		.data = NULL,
	};

	for (iout = 1; iout <= 12; iout++) {
		lsoda(&ctx, y, &t, tout, itol, rtol, atol, itask, &istate, &opt);
		printf(" at t= %12.4e y= %14.6e %14.6e %14.6e\n", t, y[1], y[2], y[3]);
		if (istate <= 0) {
			printf("error istate = %d\n", istate);
			exit(0);
		}
		tout = tout * 10.0E0;
	}

	return 0;
}
/*
 The correct answer (up to certain precision):

 at t=   4.0000e-01 y=   9.851712e-01   3.386380e-05   1.479493e-02
 at t=   4.0000e+00 y=   9.055333e-01   2.240655e-05   9.444430e-02
 at t=   4.0000e+01 y=   7.158403e-01   9.186334e-06   2.841505e-01
 at t=   4.0000e+02 y=   4.505250e-01   3.222964e-06   5.494717e-01
 at t=   4.0000e+03 y=   1.831976e-01   8.941773e-07   8.168015e-01
 at t=   4.0000e+04 y=   3.898729e-02   1.621940e-07   9.610125e-01
 at t=   4.0000e+05 y=   4.936362e-03   1.984221e-08   9.950636e-01
 at t=   4.0000e+06 y=   5.161833e-04   2.065787e-09   9.994838e-01
 at t=   4.0000e+07 y=   5.179804e-05   2.072027e-10   9.999482e-01
 at t=   4.0000e+08 y=   5.283675e-06   2.113481e-11   9.999947e-01
 at t=   4.0000e+09 y=   4.658667e-07   1.863468e-12   9.999995e-01
 at t=   4.0000e+10 y=   1.431100e-08   5.724404e-14   1.000000e+00
 */
