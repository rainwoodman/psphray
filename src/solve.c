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
	double Gamma = p[0];
	double gamma = p[1];
	double alpha = p[2];
	double nH = p[3];

	double xHI = x[0];
	double xHII = 1.0 - xHII;
	double ye = x[1];
	double yG = x[2];

/* from GABE's notes  */
	dxdt[0] = -(Gamma + gamma * ye * nH) * xHI + alpha * ye * nH * xHII;
	dxdt[1] = -dxdt[0];
	/* yG is the recombination photon */
	dxdt[2] = alpha * ye * nH * xHII;

	return GSL_SUCCESS;
}

static int jacobian(double t, const double x[], double *dfdx, double dfdt[], void * params) {
	double * p = params;
	double Gamma = p[0];
	double gamma = p[1];
	double alpha = p[2];
	double nH = p[3];

	double xHI = x[0];
	double xHII = 1.0 - xHII;
	double ye = x[1];
	double yG = x[2];

	gsl_matrix_view dfdx_mat = gsl_matrix_view_array (dfdx, 3, 3);
	gsl_matrix * m = &dfdx_mat.matrix;
	double G = (Gamma + gamma * ye * nH);
	double R = (alpha * ye * nH);
	double E = gamma * nH * xHI - alpha * nH * xHII;
	gsl_matrix_set(m, 0, 0, -G - R);
	gsl_matrix_set(m, 0, 1, -E);
	gsl_matrix_set(m, 0, 2, 0.0);
	gsl_matrix_set(m, 1, 0, G + R);
	gsl_matrix_set(m, 1, 1, E);
	gsl_matrix_set(m, 1, 2, 0.0);
	gsl_matrix_set(m, 2, 0, - alpha * ye * nH);
	gsl_matrix_set(m, 2, 1, alpha * nH * xHII);
	gsl_matrix_set(m, 2, 2, 0.0);

	dfdt[0] = 0;
	dfdt[1] = 0;
	dfdt[2] = 0;

	return GSL_SUCCESS;
}
typedef struct _Solver Solver;

struct _Solver {
	double param[4];
	gsl_odeiv2_driver * driver;
	gsl_odeiv2_system sys;
};

Solver * solver_new() {
	Solver * s = malloc(sizeof(Solver));
	s->sys.function = function;
	s->sys.jacobian = jacobian;
	s->sys.dimension = 3;
	s->sys.params = s->param;
	s->driver = gsl_odeiv2_driver_alloc_y_new(&s->sys, gsl_odeiv2_step_bsimp,
			1e-6, 1e-6, 0.0);
	return s;
}
void solver_delete(Solver * s) {
	gsl_odeiv2_driver_free(s->driver);
	free(s);
}
int solver_evolve(Solver * s, intptr_t ipar) {
	double seconds = (psys.tick - psys.lasthit[ipar]) * psys.tick_time / U_SEC;

	double logT = log10(ieye2T(psys.ie[ipar], psys.ye[ipar]));
	double Gamma = psys.deposit[ipar];
	double gamma = ar_get(AR_HI_CI, logT);
	double alpha = ar_get(AR_HII_RC_A, logT);
	double nH = C_HMF * psys.rho[ipar] / (U_MPROTON / (U_CM * U_CM * U_CM));

	s->param[0] = Gamma;
	s->param[1] = gamma;
	s->param[2] = alpha;
	s->param[3] = nH;

	double x[0];

	x[0] = psys.xHI[ipar];
	x[1] = psys.ye[ipar];
	x[2] = 0;

	double t = 0;

	gsl_odeiv2_driver_reset(s->driver);
	int code = gsl_odeiv2_driver_apply(s->driver, &t, seconds, x);

	if(code != GSL_SUCCESS) {
		ERROR("gsl failed"); /*FIXME: add a counter, and recover */
		return 0;
	}
	psys.xHI[ipar] = x[0];
	psys.ye[ipar] = x[1];
	psys.recomb[ipar] += x[2] * C_HMF * psys.mass[ipar] / (U_MPROTON);
	return 1;
}

static int test(void) {
	double Gamma = 1e-2;
	double gamma = 5.8e-11;
	double alpha = 1.3e-13;
	double nH = 10;
	double param[] = {Gamma, gamma, alpha, nH};

	gsl_odeiv2_system sys = {function, jacobian, 3, param};
	gsl_odeiv2_driver * d =
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_bsimp,
			1e-6, 1e-6, 0.0);
	double t = 0.0, t1 = 100.0;
	double x[3] = {1.0, 0.1};
	int i;
	for(i = 1; i < 1000; i++) {
		double ti = i * t1 / 100.0;
		int code;
		if(GSL_SUCCESS != (code = gsl_odeiv2_driver_apply(d, &t, ti, x))) {
			printf("gsl error, %d", code);
			exit(0);
		}
		printf("%g %g %g\n", t, x[0], x[1]);
	}
	
	gsl_odeiv2_driver_free(d);
	return 0;

}
