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
typedef struct _Solver Solver;

int solver_evolve_numerical (Solver * s, double Gamma, double gamma, double alpha, double y, double x[], double seconds);

struct _Solver {
	double param[4];
	gsl_odeiv2_driver * driver;
	gsl_odeiv2_system sys;
};

Solver * solver_new() {
	Solver * s = malloc(sizeof(Solver));
	s->sys.function = function;
	s->sys.jacobian = jacobian;
	s->sys.dimension = 1;
	s->sys.params = s->param;
	return s;
}
void solver_delete(Solver * s) {
	free(s);
}
int solver_evolve(Solver * s, intptr_t ipar) {
	if(psys.tick == psys.lasthit[ipar]) {
		ERROR("particle solved twice at one tick");
	}
	double seconds = (psys.tick - psys.lasthit[ipar]) * psys.tick_time / U_SEC;

	double logT = 4;log10(ieye2T(psys.ie[ipar], psys.ye[ipar]));
	double NH = C_HMF * psys.mass[ipar] / U_MPROTON;
	double nH = C_HMF * psys.rho[ipar] / (U_MPROTON / (U_CM * U_CM * U_CM));
	/* everything multiplied by nH, saving some calculations */
	double gamma = ar_get(AR_HI_CI, logT) * nH;
	double alpha = ar_get(AR_HII_RC_A, logT) * nH;
	double Gamma = 0;
	double y = psys.ye[ipar] - (1. - psys.xHI[ipar]);

	double x[1];
	x[0] = psys.xHI[ipar];
	int code = solver_evolve_numerical (s, Gamma, gamma, alpha, y, x, seconds);

	if(GSL_SUCCESS != code ) {
	//	WARNING("gsl failed: %s: %ld Gamma=%g gamma=%g alpha=%g nH=%g xHI=%g ye=%g, seconds=%g", 
	//		gsl_strerror(code), ipar, Gamma, gamma, alpha, nH, x[0], x[1], seconds); /*FIXME: add a counter, and recover */
		return 0;
	}

	if(x[0] > 1) x[0] = 1;
	if(x[0] < 0) x[0] = 0;

	double dxHI = x[0] - psys.xHI[ipar];
	double xHII = 1 - x[0];

	if(xHII > 1) xHII = 1;
	if(xHII < 0) xHII = 0;
	psys.xHI[ipar] = x[0];
	psys.ye[ipar] = y + xHII;
	double alpha_A = ar_get(AR_HII_RC_A, logT);

	double alpha_B = ar_get(AR_HII_RC_B, logT);
	double fac = (alpha_A - alpha_B) / alpha_B;
	psys.recomb[ipar] += dxHI * NH * fac;
	return 1;
}
int solver_evolve_analytic (Solver * s, double Gamma, double gamma, double alpha, double y, double x[], double seconds) {
	double R = (gamma + alpha);
	double Q = -(Gamma + (gamma + 2 * alpha) + (gamma + alpha) * y);
	double P = alpha * (1. + y);
	double q = - 0.5 * (Q + copysign(sqrt(Q*Q - 4 * R * P), Q));
	double x1 = P / q;
	double x2 = q / R;

	x[0] = x1;
	return GSL_SUCCESS;
}
int solver_evolve_numerical (Solver * s, double Gamma, double gamma, double alpha, double y, double x[], double seconds) {

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

	s->param[0] = R;
	s->param[1] = Q;
	s->param[2] = P;
	s->param[3] = B;

	if(seconds * d > 2e1) {
	/* if we are already equilibriem*/
	/* characteristic time is 2 / d, but the dependence is exponetial,
     * see Gabe's notes */
		x[0] = x1;
		return GSL_SUCCESS;
	}
	double t = 0;

	s->driver = gsl_odeiv2_driver_alloc_y_new(&s->sys, gsl_odeiv2_step_msbdf,
			seconds/ 10000., 1e-6, 0.0);
	gsl_odeiv2_driver_reset(s->driver);
	gsl_odeiv2_driver_set_hmin(s->driver, seconds/ 100000.);
	int code = gsl_odeiv2_driver_apply(s->driver, &t, seconds, x);

	gsl_odeiv2_driver_free(s->driver);

	if(x[0] > 1.0) {
		x[0] = 1.0;
	}
	if(x[0] < 0.0) {
		WARNING("%ld output x[0] = %e < 0.0\n", psys.tick, x[0]);
		x[0] = 0.0;
	}
	return code;
}
