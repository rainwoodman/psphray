#include <stdio.h>
#include <messages.h>
#include <stdint.h>
#include <math.h>
#include "config.h"
#include "psystem.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

int step_evolve_numerical (double Gamma, double gamma, double alpha, double y, double x[], double dx[], double seconds);

static int function(double t, const double x[], double dxdt[], void * params){
	double * p = params;
	const double R = p[0];
	const double Q = p[1];
	const double P = p[2];
	const double B = p[3];

	const double d = x[0] - B;
/* from GABE's notes  */
	dxdt[0] = R * x[0] * x[0] + Q * x[0] + P;
	dxdt[0] *= - tanh(10 * d);
/* force the unstable solution to converge back to the stable one */

	return GSL_SUCCESS;
}

static int jacobian(double t, const double x[], double *dfdx, double dfdt[], void * params) {
	double * p = params;
	const double R = p[0];
	const double Q = p[1];
	const double P = p[2];
	const double B = p[3];

	const double d = x[0] - B;

	dfdx[0] = 2 * R * x[0] + Q;
/* force the unstable solution to converge back to the stable one */
	const double th = tanh(- 10 * d);
	dfdx[0] += - (R * x[0] * x[0] + Q * x[0] + P) *  (1 - th * th) * 10;
	dfdt[0] = 0;

	return GSL_SUCCESS;
}


int step_evolve(Step * step, double time) {
	const double seconds = time / U_SEC;

	const double T = step->T;
	const double logT = log10(T);
	const double gamma_HI = ar_get(AR_HI_CI, logT) * step->nH;
	const double alpha_HII = ar_get(AR_HII_RC_A, logT) * step->nH;
	const double alpha_HII1 = alpha_HII - ar_get(AR_HII_RC_B, logT) * step->nH;
	const double Gamma_HI = 0;

	double x[1];
	x[0] = step->xHI;
	const double y = step->y;
	
	double dx[1] = {0};

	const int code = step_evolve_numerical (Gamma_HI, gamma_HI, alpha_HII, y, x, dx, seconds);

	if(GSL_SUCCESS != code ) {
		WARNING("gsl failed: %s: T=%g, Gamma=%g gamma=%g alpha=%g y=%g xHI=%g, seconds=%g", 
			gsl_strerror(code), T, Gamma_HI, gamma_HI, alpha_HII,  y, x[0], seconds);
		return 0;
	}

	const double fac = (alpha_HII1 / alpha_HII);
	step->dxHI = dx[0];
	step->dye = - dx[0];
	step->dyGH = dx[0] * fac;

	return 1;
}
int step_evolve_analytic (double Gamma, double gamma, double alpha, double y, double x[], double dx[], double seconds) {
	const double R = (gamma + alpha);
	const double Q = -(Gamma + (gamma + 2 * alpha) + (gamma + alpha) * y);
	const double P = alpha * (1. + y);
	const double q = - 0.5 * (Q + copysign(sqrt(Q*Q - 4 * R * P), Q));
	if(q == 0.0) {
		/* assuming this is due to gamma -> 0 slower than alpha->0.
		 * in that case R and P also -> 0 and a good approximation is 1 */
		dx[0] = 1.0 - x[0];
		return GSL_SUCCESS;
	}
	const double x1 = P / q;
	const double x2 = q / R;

	dx[0] = x1 - x[0];
	return GSL_SUCCESS;
}
int step_evolve_numerical (double Gamma, double gamma, double alpha, double y, double x[], double dx[], double seconds) {

	const double R = (gamma + alpha);
	const double Q = -(Gamma + (gamma + 2 * alpha) + (gamma + alpha) * y);
	const double P = alpha * (1. + y);
	const double d = sqrt(Q*Q - 4 * R * P);
	const double q = - 0.5 * (Q + copysign(d, Q));
	const double x1 = P / q;
	const double x2 = q / R;
	const double B = (x2 > x1)? x2: x1;

//	MESSAGE("x1 = %le , x2 = %le x[0] - B %g", x1, x2, x[0] - B);

	gsl_odeiv2_system sys;
	gsl_odeiv2_driver * driver;

	double param[4] = {R * seconds, Q * seconds, P * seconds, B};
	sys.function = function;
	sys.jacobian = jacobian;
	sys.dimension = 1;
	sys.params = param;

	if(d > 2e1) {
	/* if we are already equilibriem*/
	/* characteristic time is 2 / d, but the dependence is exponetial,
     * see Gabe's notes */
		dx[0] = x1 - x[0];
		return GSL_SUCCESS;
	}
	/* when the character time is too big, assume nothing happened.*/
/*
	if(gamma * seconds < 1e-2 && alpha * seconds < 1e-2) {
		dx[0] = 0;
		return GSL_SUCCESS;
	}
*/

	double t = 0;


	dx[0] = x[0];
	driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msbdf,
			1.0 / 10000., 1e-6, 0.0);
//	gsl_odeiv2_driver_reset(driver);
	//gsl_odeiv2_driver_set_hmin(driver, seconds/ 100000.);
	int code = gsl_odeiv2_driver_apply(driver, &t, 1.0, dx);

	gsl_odeiv2_driver_free(driver);

	if(dx[0] > 1.0) {
		dx[0] = 1.0;
	}
	if(dx[0] < 0.0) {
		WARNING("%ld output x[0] = %e < 0.0\n", psys.tick, dx[0]);
		dx[0] = 0.0;
	}
	dx[0] -= x[0];
	return code;
}
