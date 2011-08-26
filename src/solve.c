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
	const double Y = p[3];
	const double Z = p[4];

	double xHII, xHI;
	lambdaHI_to_xHI_xHII(x[0], xHI, xHII);

/* from GABE's notes  */
	const double fac = xHI * xHI + xHII * xHII;
	dxdt[0] = (R * xHI * xHI + Q * xHI + P) / fac;
	dxdt[1] = (Z * xHII * xHII + Y * xHII);
	if(xHI >= 1.0) {
		dxdt[0] *= -1;
		dxdt[1] *= -1;
	}
	return GSL_SUCCESS;
}

static int jacobian(double t, const double x[], double *dfdx, double dfdt[], void * params) {
	double * p = params;
	const double R = p[0];
	const double Q = p[1];
	const double P = p[2];
	const double Y = p[3];
	const double Z = p[4];

	double xHI, xHII;
	lambdaHI_to_xHI_xHII(x[0], xHI, xHII);

	const double fac = xHI * xHI + xHII * xHII;

	dfdx[0] = 2 * R * xHI + Q + (R * xHI * xHI + Q * xHI + P) * (2 * xHII - 2 * xHI) / fac;

	dfdx[1] = 0;
	dfdx[2] = - fac * (2 * Z * xHII + Y);
	dfdx[3] = 0;
	dfdt[0] = 0;
	dfdt[1] = 0;

	if(xHI >= 1.0) {
		dfdx[0] *= -1;
		dfdx[2] *= -1;
	}
	return GSL_SUCCESS;
}


int step_evolve(Step * step, double time) {
	const double seconds = time / U_SEC;

	const double T = step->T;
	const double logT = log10(T);
	const double gamma_HI = ar_get(AR_HI_CI, logT) * step->nH;
	const double alpha_HII = ar_get(AR_HII_RC_A, logT) * step->nH;
	const double alpha_HII1 = alpha_HII - ar_get(AR_HII_RC_B, logT) * step->nH;
	if(alpha_HII1 < 0) abort();
	const double Gamma_HI = 0;

	double x[2];
	x[0] = step->lambdaHI;
	x[1] = 0;
	const double yeMET = step->yeMET;
	
	double dx[1] = {0};

	const double R = (gamma_HI + alpha_HII);
	const double Q = -(Gamma_HI + (gamma_HI + 2 * alpha_HII) + (gamma_HI + alpha_HII) * yeMET);
	const double P = alpha_HII * (1. + yeMET);
	const double Y = alpha_HII1 * yeMET;
	const double Z = alpha_HII1;

//	MESSAGE("x1 = %le , x2 = %le x[0] - B %g", x1, x2, x[0] - B);

	gsl_odeiv2_system sys;
	gsl_odeiv2_driver * driver;

	double param[5] = {R * seconds, Q * seconds, P * seconds, Y * seconds, Z * seconds};
	sys.function = function;
	sys.jacobian = jacobian;
	sys.dimension = 2;
	sys.params = param;

	double t = 0;

	driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msbdf,
			1.0 / 10000., 1e-6, 1e-6);
//	gsl_odeiv2_driver_reset(driver);
	//gsl_odeiv2_driver_set_hmin(driver, seconds/ 100000.);
	int code = gsl_odeiv2_driver_apply(driver, &t, 1.0, x);

	gsl_odeiv2_driver_free(driver);

	if(GSL_SUCCESS != code ) {
		WARNING("gsl failed: %s: T=%g, Gamma=%g gamma=%g alpha=%g yeMET =%g xHI=%g, seconds=%g", 
			gsl_strerror(code), T, Gamma_HI, gamma_HI, alpha_HII,  yeMET, x[0], seconds);
		return 0;
	}

	if(x[0] < 0) x[0] = 0;

	step->dyGrec = x[1];
	step->lambdaHI = x[0];

	return 1;
}

int step_evolve_numerical (double Gamma, double gamma, double alpha, double y, double x[], double dx[], double seconds) {

}
