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

/* ******************
 * Equations:
 * ==============
 *
 * for x[0]:
 * d lambda   1
 * ------ = ---- R x_I ** 2 + Q * x_I + P  [fac = x_I **2 + x_II **2]
 * d t       fac
 *
 * for x[1]:
 * d yGrec
 * ------- = Z x_II ** 2 + Z * x_II
 * d t
 *
 * for x[2]:
 * d u         1
 * ---   =    --- (H - L)
 * d t        rho
 */

/*  
 * factor all coefficients by corresponding (nH) powers
 * ==============
 * gamma_HI alpha_HII shall take (nH) ** 2, but one of the powers
 * is cancelled with the nH on the left-hand-side of the equation
 *
 * eta is nH ** 2, psi_HI is nH **2, psi_HeI is nH ** 3 etc 
 * */
#define FETCH_VARS \
	const double seconds = step->time / U_SEC; \
	const double T = step->T; \
	const double logT = log10(T); \
	const double gamma = seconds * ar_get(AR_HI_CI, logT) * step->nH; \
	const double alpha_A = seconds * ar_get(AR_HII_RC_A, logT) * step->nH; \
	const double alpha_AB = alpha_A - seconds * ar_get(AR_HII_RC_B, logT) * step->nH; \
	const double nH2 = step->nH * step->nH; \
	const double nH3 = nH2 * step->nH; \
	const double eta_HII = seconds * ar_get(AR_HII_RCC_A, logT) * nH2; \
	const double psi_HI = seconds * ar_get(AR_HI_CEC_A, logT) * nH2; \
	const double zeta_HI = seconds * ar_get(AR_HI_CIC_A, logT) * nH2; \
\
	const double yGdep_mean = step->yGdep; \
	const double yeMET = step->yeMET; \
\
	double xHII, xHI; \
	lambdaHI_to_xHI_xHII(x[0], xHI, xHII); \
	const double fac = xHI * xHI + xHII * xHII; \
	const double ye = xHII + step->yeMET; \
	const int sign = xHI >= 1.0? -1: 1;

static int function(double t, const double x[], double dxdt[], Step * step){

	FETCH_VARS;

/* from GABE's notes  */
	dxdt[0] = sign * (-yGdep_mean - gamma * ye * xHI + alpha_A* ye * xHII) / fac;
	dxdt[1] = sign * alpha_AB * ye * xHII;
	dxdt[2] = 0.0;

	return GSL_SUCCESS;
}

static int jacobian(double t, const double x[], double *dfdx, double dfdt[], Step * step) {

	const int D = 3;
	FETCH_VARS;

	const double dxdt0 = (-yGdep_mean - gamma * ye * xHI + alpha_A * ye * xHII) / fac;
	dfdx[0 * D + 0] = sign * (-2 * (xHI - xHII) * dxdt0 + gamma * xHI - (gamma * ye + alpha_A * ye + alpha_A * xHII));

	dfdx[0 * D + 1] = 0;
	dfdx[0 * D + 2] = 0;
	dfdx[1 * D + 0] = sign * (-fac * alpha_AB * (ye + xHII));
	dfdx[1 * D + 1] = 0;
	dfdx[1 * D + 2] = 0;
	dfdx[2 * D + 0] = 0;
	dfdx[2 * D + 1] = 0;
	dfdx[2 * D + 2] = 0;

	dfdt[0] = 0;
	dfdt[1] = 0;
	dfdt[2] = 0;
	return GSL_SUCCESS;
}


int step_evolve(Step * step) {

//	MESSAGE("x1 = %le , x2 = %le x[0] - B %g", x1, x2, x[0] - B);

	gsl_odeiv2_system sys;
	gsl_odeiv2_driver * driver;

	sys.function = function;
	sys.jacobian = jacobian;
	sys.dimension = 3;
	sys.params = step;

	double t = 0;

	double x[3];
	x[0] = step->lambdaHI;
	x[1] = 0;
	x[2] = step->ie;

	driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msbdf,
			1.0 / 100., 1e-6, 1e-6);
//	gsl_odeiv2_driver_reset(driver);
	gsl_odeiv2_driver_set_nmax(driver, 100000);
	int code = gsl_odeiv2_driver_apply(driver, &t, 1.0, x);

	gsl_odeiv2_driver_free(driver);

	if(GSL_SUCCESS != code ) {
		WARNING("gsl failed: %s", gsl_strerror(code));
		return 0;
	}

	if(x[0] < 0) x[0] = 0;

	step->lambdaHI = x[0];
	step->dyGrec = x[1];
	step->ie = x[2];

	return 1;
}

int step_evolve_numerical (double Gamma, double gamma, double alpha, double y, double x[], double dx[], double seconds) {

}
