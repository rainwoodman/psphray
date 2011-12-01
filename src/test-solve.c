#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <messages.h>
#include "reader.h"
#include "config.h"
#include "psystem.h"


extern PSystem psys;
typedef struct _Solver Solver;
extern int step_solve(Step * step);

int main(int argc, char* argv[]) {
	cfg_init("../test/configfile");

	CFG_ISOTHERMAL=1;
	CFG_ADIABATIC=0;
	CFG_ON_THE_SPOT=0;
	CFG_SECONDARY_IONIZATION = 0;
	CFG_H_ONLY = 1;
	Step step = {0};
step.T = 1e4;
step.lambdaH = 1.0;
step.nH = 1.786e67;
step.yeMET= 0;
step.ie = 1e4;
step.yGdepHI = 0.000001;

	double RC, CI, RCB;
	RC=ar_get(AR_HII_RC_A, log10(step.T));
	RCB = ar_get(AR_HII_RC_B, log10(step.T));
	CI=ar_get(AR_HI_CI, log10(step.T));
	printf("1e4K, ie= %g, check=%g\n", Tye2ie(1e4, 0), ieye2T(Tye2ie(1e4, 0), 0));
	printf("%g K, RC_A = %g RC_B = %g, CI = %g\n", 
		step.T,
		RC / (U_CM * U_CM * U_CM / U_SEC),
		RCB / (U_CM * U_CM * U_CM / U_SEC),
		CI / (U_CM * U_CM * U_CM / U_SEC)
	);
	double P, Q, R;
	P = RC;
	Q = -(CI + 2 * RC);
	R = CI + RC;
	double q = -0.5 * (Q + copysign(sqrt(Q*Q - 4 * R*P), Q));
	
	printf("analytical limite=%g\n", P / q);
		

	printf("recomb time = %g myr\n", 1.0 / (step.nH * RCB) / U_MYR);
	double time = 0.01 * U_MYR;
	intptr_t i;
	step.time = time;
	Step s2 = step;
	for(i = 0; i < 1001; i++) {
		const double xHI = lambdaH_to_xHI(step.lambdaH);
		const double xHII = lambdaH_to_xHII(step.lambdaH);
		const double xHI2 = lambdaH_to_xHI(s2.lambdaH);
		const double xHII2 = lambdaH_to_xHII(s2.lambdaH);
		printf("%g %g %g %g %g %g %g\n", i * time / U_MYR,  xHI, xHI2, xHII, xHII2, step.yGrecHII, s2.yGrecHII);

		int code = step_evolve_euler (&step);
		if(!code) {
			WARNING("evolve failed");
		}
		code = step_evolve_lsoda (&s2);
		if(!code) {
			WARNING("evolve failed");
		}
	}

	return 0;
}
