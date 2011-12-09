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
	const double xHI = lambdaH_to_xHI(x[0]); \
	const double xHII = lambdaH_to_xHII(x[0]); \
	const double xHeI = lambdaHe_to_xHeI(x[1], x[1]); \
	const double xHeII = lambdaHe_to_xHeII(x[1], x[2]); \
	const double xHeIII = lambdaHe_to_xHeIII(x[1], x[2]); \
	const double fac = xHI * xHI + xHII * xHII; \
	const double ye = lambda_to_ye(x[0], x[1], x[2]) + step->yeMET; \
	const double nH2 = step->nH * step->nH; \
	const double nH3 = nH2 * step->nH; \
	const double nH = step->nH; \
	const double yeMET = step->yeMET; \
	const double T = CFG_ISOTHERMAL?step->T:ieye2T(x[6], ye); \
	const double logT = log10(T); \
	const double gamma_HI = ar_get(AR_HI_CI, logT) * nH; \
	const double alpha_HII_A = ar_get(AR_HII_RC_A, logT) * nH; \
	const double alpha_HII_B = ar_get(AR_HII_RC_B, logT) * nH; \
	const double alpha_HII_AB = fdim(alpha_HII_A, alpha_HII_B); \
	const double gamma_HeI = ar_get(AR_HEI_CI, logT) * nH; \
	const double gamma_HeII = ar_get(AR_HEII_CI, logT) * nH; \
	const double alpha_HeII_A = ar_get(AR_HEII_RC_A, logT) * nH; \
	const double alpha_HeII_B = ar_get(AR_HEII_RC_B, logT) * nH; \
	const double alpha_HeII_AB = fdim(alpha_HeII_A, alpha_HeII_B); \
	const double alpha_HeIII_A = ar_get(AR_HEIII_RC_A, logT) * nH; \
	const double alpha_HeIII_B = ar_get(AR_HEIII_RC_B, logT) * nH; \
	const double alpha_HeIII_AB = fdim(alpha_HeIII_A, alpha_HeIII_B); \
	const double eta_HII = (CFG_ON_THE_SPOT?ar_get(AR_HII_RCC_B, logT):ar_get(AR_HII_RCC_A, logT)) * nH; \
	const double psi_HI = ar_get(AR_HI_CEC, logT) * nH; \
	const double zeta_HI = ar_get(AR_HI_CIC, logT) * nH; \
	const double beta = ar_get(AR_E_BREMC, logT) * nH; \
	const double chi = ar_get(AR_E_COMPC, logT) ; \
\
\

static int function(double t, const double x[], double dxdt[], Step * step){

	int i;

	FETCH_VARS;

/* from GABE's notes  */
	if(!CFG_ON_THE_SPOT) {
		dxdt[0] = -step->yGdepHI / step->time + (- gamma_HI * ye * xHI + alpha_HII_A * ye * xHII);
		dxdt[3] = alpha_HII_AB * ye * xHII;
		if(CFG_H_ONLY) {
			dxdt[1] = 0.0;
			dxdt[2] = 0.0;
			dxdt[4] = 0.0;
			dxdt[5] = 0.0;
		} else {
			dxdt[1] = -step->yGdepHeI / step->time + (- gamma_HeI * ye * xHeI + alpha_HeII_A * ye * xHeII);
			dxdt[2] = -step->yGdepHeII / step->time + (- gamma_HeII * ye * xHeII - alpha_HeII_A * ye * xHeII + alpha_HeIII_A * ye * xHeIII);
			dxdt[4] = alpha_HeII_AB * ye * xHeII;
			dxdt[5] = alpha_HeIII_AB * ye * xHeIII;
		}
	} else {
		dxdt[0] = -step->yGdepHI / step->time + (- gamma_HI * ye * xHI + alpha_HII_B * ye * xHII);
		if(CFG_H_ONLY) {
			dxdt[1] = 0.0;
			dxdt[2] = 0.0;
		} else {
			dxdt[1] = -step->yGdepHeI / step->time + (-gamma_HeI * ye * xHI + alpha_HeII_B * ye * xHeII);
			dxdt[2] = -step->yGdepHeII / step->time + (-gamma_HeII * ye * xHI - alpha_HeII_B * ye * xHeII + alpha_HeIII_B * ye * xHeIII);
		}
		dxdt[3] = 0;
		dxdt[4] = 0;
		dxdt[5] = 0;
	}
	if(CFG_ADIABATIC | CFG_ISOTHERMAL) {
		dxdt[6] = 0.0;
	} else {
		const double L = C_H_PER_MASS * (
			(zeta_HI + psi_HI) * ye * xHI + 
			 (eta_HII * ye * xHII) +
			 beta * ye * xHII + chi * ye);
	//	MESSAGE("T = %g H=%g L=%g eta=%g psi=%g zeta=%g beta=%g xHI=%g\n", T, heat_mean, L, eta_HII, psi_HI, zeta_HI, beta, xHI);
		dxdt[6] = step->heat / step->time - L;
	}
	for(i = 0; i < 7; i++) {
		dxdt[i] *= step->time;
		if(isinf(dxdt[i]) || isnan(dxdt[i])) {
			ERROR("nan found in dxdt[%d]", i);
		}
	}
/*
	if(dxdt[3] < 0) {
		ERROR("dxdt[3] < 0");
	}
	if(dxdt[4] < 0) {
		ERROR("dxdt[4] < 0");
	}
	if(dxdt[5] < 0) {
		ERROR("dxdt[5] < 0");
	}
*/
	if(x[6] + dxdt[6] < 0) {
		ERROR("x[6] < 0");
	}
	return GSL_SUCCESS;
}

#include <lsoda/lsoda.h>
struct step_evolve_t {
	struct lsoda_context_t ctx;
	struct lsoda_opt_t opt;
};
static int mxstep = 1000;
void * step_evolve_prepare() {
	struct step_evolve_t * control;
	static double rtol[7] = {1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4};
	static double atol[7] = {1e-7, 1e-7, 1e-7, 1e-7, 1e-7, 1e-7, 1.};

	control = calloc(1, sizeof(struct step_evolve_t));
	control->ctx.function = function;
	control->ctx.neq = 7;
	control->opt.ixpr = 0;
	control->opt.mxstep = mxstep;
	control->opt.rtol = rtol;
	control->opt.atol = atol;
	control->opt.itask = 1;
	lsoda_prepare(&control->ctx, &control->opt);
	return control;
}
void step_evolve_free(void * control) {
	struct step_evolve_t * e = control;
	lsoda_free(&e->ctx);
	free(control);
}
int step_evolve_lsoda(void * control, Step * step) {
	struct step_evolve_t * e = control;
	e->ctx.data = step;
	e->ctx.state = 1;

	lsoda_reset(&e->ctx);

	double x[7];
	double dxdt[7];
	double t = 0;
	x[0] = step->lambdaH;
	x[1] = step->lambdaHeI;
	x[2] = step->lambdaHeII;
	x[3] = step->yGrecHII;
	x[4] = step->yGrecHeII;
	x[5] = step->yGrecHeIII;
	x[6] = step->ie;
	do {
		e->ctx.state = 1;
		lsoda(&e->ctx, x, &t, 1.);
		if(e->ctx.error) {
			WARNING("lsodaerror %s\n", e->ctx.error);
			free(e->ctx.error);
			e->ctx.error = NULL;
		}
		if(e->ctx.state == -1) {
			FETCH_VARS;
			e->opt.mxstep = e->opt.mxstep * 2;
			if(mxstep < e->opt.mxstep) mxstep = e->opt.mxstep;
			WARNING("continuing with new mxstep = %d", e->opt.mxstep);
			MESSAGE("%lu %ld %g %g %g %g %g %g %g %g %g %g\n",
						psys.tick, step->ipar, T, xHI, xHII, xHeI, xHeII, xHeIII, gamma_HI, alpha_HII_A, alpha_HeII_A, alpha_HeIII_A);
		}
	} while(e->ctx.state == -1 && t < 1.);

	if(e->ctx.state != 2) {
		ERROR("lsoda failed");
	}
	x[0] = fmin(1., fmax(0., x[0]));
	x[1] = fmin(1., fmax(0., x[1]));
	x[2] = fmin(1. - x[1], fmax(0., x[2]));
	x[3] = fmax(0., x[3]);
	x[4] = fmax(0., x[4]);
	x[5] = fmax(0., x[5]);

	step->lambdaH = x[0];
	step->lambdaHeI = x[1];
	step->lambdaHeII = x[2];
	step->yGrecHII = x[3];
	step->yGrecHeII = x[4];
	step->yGrecHeIII = x[5];
	step->ie = x[6];
	step->step_remain = 0.0;
	return 1;
}
int step_evolve_euler(void * control, Step * step) {

//	MESSAGE("x1 = %le , x2 = %le x[0] - B %g", x1, x2, x[0] - B);

	double t = 0;

	double x[7];
	double dxdt[7];
	x[0] = step->lambdaH;
	x[1] = step->lambdaHeI;
	x[2] = step->lambdaHeII;
	x[3] = step->yGrecHII;
	x[4] = step->yGrecHeII;
	x[5] = step->yGrecHeIII;
	x[6] = step->ie;

	if(x[3] > 1.0) {
		ERROR("ygrecHII too big to start with");
	}
	FETCH_VARS;

	double rectime = 1.0 / (alpha_HII_A);
	double iontime = 1.0 / (gamma_HI);
	if(!CFG_H_ONLY) {
		rectime = fmin(rectime, 1.0 / alpha_HeII_A);
		rectime = fmin(rectime, 1.0 / alpha_HeIII_A);
		iontime = fmin(iontime, 1.0 / gamma_HeI);
		iontime = fmin(iontime, 1.0 / gamma_HeII);
	}
	double internal_time = fmin(rectime, iontime) / 200.;
	internal_time = fmin(internal_time, step->time);

	int nsteps = step->time / internal_time;
	if(nsteps > 10000) nsteps = 10000;
	double internal_step = 1.0 / nsteps;

	if(nsteps > 1) step->refined = 1;
	int k;
	for(k = 0; k < nsteps; k ++) {
		t = k * internal_step;
		int i;
		for(i = 0 ; i < 7; i++) {
			if(isinf(x[i]) || isnan(x[i])) {
				ERROR("nan found in x[%d]", i);
			}
		}

		function(t, x, dxdt, step);

		x[0] += dxdt[0] * internal_step;
		x[1] += dxdt[1] * internal_step;
		x[2] += dxdt[2] * internal_step;
		x[3] += dxdt[3] * internal_step;
		x[4] += dxdt[4] * internal_step;
		x[5] += dxdt[5] * internal_step;
		x[6] += dxdt[6] * internal_step;
		for(i = 0 ; i < 7; i++) {
			if(isinf(x[i]) || isnan(x[i])) {
				ERROR("nan2 found in x[%d]", i);
			}
		}
		for(i = 0 ; i < 3; i++) {
			if(x[i] < 0.0 || x[i] > 1.0) {
				k++;
				goto terminate;
			}
		}
		/* need to check xHeII too*/
		if(1.0 - x[1] - x[2] < 0.0) {
			k++;
			goto terminate;
		}
	}
	terminate:
//	for(i = 0; i <= 5; i++) {
//		x[i] = fmax(0, fmin(1, x[i]));
//	}
	if(x[3] > 1.0) {
		ERROR("yGrecHII grows > 1.0 after ");
	}

	step->lambdaH = x[0];
	step->lambdaHeI = x[1];
	step->lambdaHeII = x[2];
	step->yGrecHII = x[3];
	step->yGrecHeII = x[4];
	step->yGrecHeIII = x[5];
	step->ie = x[6];

	step->step_remain = 1.0 - (float) k / nsteps;
	return 1;
}
int step_evolve(void * control, Step * step) {
#ifdef USE_LSODA
	return step_evolve_lsoda(control, step);
#else
	return step_evolve_euler(control, step);
#endif
}
