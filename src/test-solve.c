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
	CFG_ON_THE_SPOT=1;
	double RC, CI;
	printf("1e4K, ie= %g, check=%g\n", Tye2ie(1e4, 0), ieye2T(Tye2ie(1e4, 0), 0));
	printf("1e4K, RC_A = %g RC_B = %g, CI = %g\n", 
		RC=ar_get(AR_HII_RC_A, log10(1e4)), 
		ar_get(AR_HII_RC_B, log10(1e4)),
		CI=ar_get(AR_HI_CI, log10(1e4))
	);
	double P, Q, R;
	P = RC;
	Q = -(CI + 2 * RC);
	R = CI + RC;
	double q = -0.5 * (Q + copysign(sqrt(Q*Q - 4 * R*P), Q));
	
	printf("analytical limite=%g\n", P / q);
		
	Step step = {0};

step.T = 1e4 + 100;
step.lambdaH = 1.569;
step.nH = 1e-3;
step.yeMET=0;
step.ie = 123;

	double time = 5000.0 * U_MYR;
	intptr_t i;
	for(i = 0; i < 10; i++) {
		const double xHI = lambdaH_to_xHI(step.lambdaH);
		const double xHII = lambdaH_to_xHII(step.lambdaH);
		printf("%g %g %g %g %g %g \n", i * time,  step.ie, step.lambdaH, xHI, xHII, step.yGrecHII);

		step.heat = 0.0;
		step.yGrecHII = 0.0;
		step.time = time;
		int code = step_evolve (&step);
		if(!code) {
			WARNING("evolve failed");
		}
	}

	return 0;
}
