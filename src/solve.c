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
	double xHII, xHI; \
	lambdaHI_to_xHI_xHII(x[0], xHI, xHII); \
	const double fac = xHI * xHI + xHII * xHII; \
	const double ye = xHII + step->yeMET; \
	const double nH2 = step->nH * step->nH; \
	const double nH3 = nH2 * step->nH; \
	const double nH = step->nH; \
	const double yGdep_mean = step->yGdep / seconds; \
	const double heat_mean = step->heat / seconds; \
	const double yeMET = step->yeMET; \
	const double T = CFG_ISOTHERMAL?step->T:ieye2T(x[2], ye); \
	const double logT = log10(T); \
	const double gamma = ar_get(AR_HI_CI, logT) * nH; \
	const double alpha_A = ar_get(AR_HII_RC_A, logT) * nH; \
	const double alpha_B = ar_get(AR_HII_RC_B, logT) * nH; \
	const double alpha_AB = alpha_A - alpha_B; \
	const double eta_HII = (CFG_ON_THE_SPOT?ar_get(AR_HII_RCC_B, logT):ar_get(AR_HII_RCC_A, logT)) * nH; \
	const double psi_HI = ar_get(AR_HI_CEC, logT) * nH; \
	const double zeta_HI = ar_get(AR_HI_CIC, logT) * nH; \
	const double beta = ar_get(AR_E_BREMC, logT) * nH; \
	const double chi = ar_get(AR_E_COMPC, logT) ; \
\
\
	const int sign = xHI >= 1.0? -1: 1;

static int function(double t, const double x[], double dxdt[], Step * step){

	FETCH_VARS;

/* from GABE's notes  */
	if(!CFG_ON_THE_SPOT) {
		dxdt[1] = sign * seconds * alpha_AB * ye * xHII;
		dxdt[0] = sign * seconds * (-yGdep_mean - gamma * ye * xHI + alpha_A * ye * xHII) / fac;
	} else {
		dxdt[1] = 0;
		dxdt[0] = sign * seconds * (-yGdep_mean - gamma * ye * xHI + alpha_B * ye * xHII) / fac;
	}
	if(CFG_ADIABATIC | CFG_ISOTHERMAL) {
		dxdt[2] = 0.0;
	} else {
		const double L = U_ERG * C_H_PER_MASS * (
			(zeta_HI + psi_HI) * ye * xHI + 
			 (eta_HII * ye * xHII) +
			 beta * ye * xHII + chi * ye);
	//	MESSAGE("T = %g H=%g L=%g eta=%g psi=%g zeta=%g beta=%g xHI=%g\n", T, heat_mean, L, eta_HII, psi_HI, zeta_HI, beta, xHI);
		dxdt[2] = sign * seconds * (heat_mean   - L);
	}

	return GSL_SUCCESS;
}

static int jacobian(double t, const double x[], double *dfdx, double dfdt[], Step * step) {

	const int D = 3;
	FETCH_VARS;

	if(!CFG_ON_THE_SPOT) {
		const double dxdt0 = (-yGdep_mean - gamma * ye * xHI + alpha_A * ye * xHII) / fac;
		dfdx[0 * D + 0] = sign * seconds * (
			- 2 * (xHI - xHII) * dxdt0 
			+ gamma * xHI 
			- ((gamma + alpha_A) * ye + alpha_A * xHII)
			);
	} else {
		const double dxdt0 = (-yGdep_mean - gamma * ye * xHI + alpha_B * ye * xHII) / fac;
		dfdx[0 * D + 0] = sign * seconds * (
			- 2 * (xHI - xHII) * dxdt0 
			+ gamma * xHI 
			- ((gamma + alpha_B) * ye + alpha_B * xHII)
			);

	}
	dfdx[0 * D + 1] = 0;
	dfdx[0 * D + 2] = 0;
	if(!CFG_ON_THE_SPOT) {
		dfdx[1 * D + 0] = 0;
	} else {
		dfdx[1 * D + 0] = sign * seconds * (-fac * alpha_AB * (ye + xHII));
	}
	dfdx[1 * D + 1] = 0;
	dfdx[1 * D + 2] = 0;
	if(CFG_ADIABATIC | CFG_ISOTHERMAL) {
		dfdx[2 * D + 0] = 0;
	} else {
		dfdx[2 * D + 0] = sign * seconds * U_ERG * C_H_PER_MASS * fac * (
			+ (eta_HII + beta) * (ye  + xHII)
			+ (zeta_HI + psi_HI) * (xHI - ye)
			+ chi
			);
	}
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
			1.0 / 100., 1e-7, 1e-7);
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
