#include <stdio.h>
#include <messages.h>
#include <stdint.h>
#include <math.h>
#include "config.h"
#include "psystem.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

extern PSystem psys;

static int function(double t, const double x[], double dxdt[], void * params){
	double * p = params;
	double R = p[0];
	double Q = p[1];
	double P = p[2];
	double B = p[3];

	double d = x[0] - B;
/* from GABE's notes  */
	dxdt[0] = R * x[0] * x[0] + Q * x[0] + P;
	dxdt[0] *= - tanh(10 * d);
/* force the unstable solution to converge back to the stable one */

	return GSL_SUCCESS;
}

static int jacobian(double t, const double x[], double *dfdx, double dfdt[], void * params) {
	double * p = params;
	double R = p[0];
	double Q = p[1];
	double P = p[2];
	double B = p[3];

	double d = x[0] - B;

	dfdx[0] = 2 * R * x[0] + Q;
/* force the unstable solution to converge back to the stable one */
	double th = tanh(- 10 * d);
	dfdx[0] += - (R * x[0] * x[0] + Q * x[0] + P) *  (1 - th * th) * 10;
	dfdt[0] = 0;

	return GSL_SUCCESS;
}


int step_evolve(Step * step, double time) {
	double seconds = time / U_SEC;

	double logT = log10(step->T);
	double gamma_HI = ar_get(AR_HI_CI, logT) * step->nH;
	double alpha_HII = ar_get(AR_HII_RC_A, logT) * step->nH;
	double alpha_HII1 = alpha_HII - ar_get(AR_HII_RC_B, logT) * step->nH;
	double Gamma_HI = 0;

	double x[1];
	x[0] = step->xHI;
	double y = step->y;

	double dx[1] = {0};
	
	int code = step_evolve_numerical (Gamma_HI, gamma_HI, alpha_HII, y, x, dx, seconds);

	if(GSL_SUCCESS != code ) {
	//	WARNING("gsl failed: %s: %ld Gamma=%g gamma=%g alpha=%g nH=%g xHI=%g ye=%g, seconds=%g", 
	//		gsl_strerror(code), ipar, Gamma, gamma, alpha, nH, x[0], x[1], seconds); /*FIXME: add a counter, and recover */
		return 0;
	}

	double fac = (alpha_HII1 / alpha_HII);

	step->dxHI = dx[0];
	step->dye = - dx[0];
	step->dyGH = dx[0] * fac;

	return 1;
}
int step_evolve_analytic (double Gamma, double gamma, double alpha, double y, double x[], double dx[], double seconds) {
	double R = (gamma + alpha);
	double Q = -(Gamma + (gamma + 2 * alpha) + (gamma + alpha) * y);
	double P = alpha * (1. + y);
	double q = - 0.5 * (Q + copysign(sqrt(Q*Q - 4 * R * P), Q));
	double x1 = P / q;
	double x2 = q / R;

	dx[0] = x1 - x[0];
	return GSL_SUCCESS;
}
int step_evolve_numerical (double Gamma, double gamma, double alpha, double y, double x[], double dx[], double seconds) {

	double R = (gamma + alpha);
	double Q = -(Gamma + (gamma + 2 * alpha) + (gamma + alpha) * y);
	double P = alpha * (1. + y);
	double d = sqrt(Q*Q - 4 * R * P);
	double q = - 0.5 * (Q + copysign(d, Q));
	double x1 = P / q;
	double x2 = q / R;
	double B;
	if(x2 > x1) B = x2;
	else B = x1;
//	MESSAGE("x1 = %le , x2 = %le x[0] - B %g", x1, x2, x[0] - B);

	gsl_odeiv2_system sys;
	gsl_odeiv2_driver * driver;

	double param[4] = {R, Q, P, B};
	sys.function = function;
	sys.jacobian = jacobian;
	sys.dimension = 1;
	sys.params = param;

	if(seconds * d > 2e1) {
	/* if we are already equilibriem*/
	/* characteristic time is 2 / d, but the dependence is exponetial,
     * see Gabe's notes */
		dx[0] = x1 - x[0];
		return GSL_SUCCESS;
	}
	double t = 0;


	dx[0] = x[0];
	driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msbdf,
			seconds/ 10000., 1e-6, 0.0);
	gsl_odeiv2_driver_reset(driver);
	gsl_odeiv2_driver_set_hmin(driver, seconds/ 100000.);
	int code = gsl_odeiv2_driver_apply(driver, &t, seconds, dx);

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
