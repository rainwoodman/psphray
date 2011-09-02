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

	printf("1e4K, ie= %g, check=%g\n", Tye2ie(1e4, 0), ieye2T(Tye2ie(1e4, 0), 0));
	Step step = {0};
step.T = 99.8801; 
step.lambdaHI=1.56959; 
step.yeMET=0;
step.nH=0.00105194;
step.yGdep=3.14168e-05;
step.ie=1.23988;
step.heat=0.0949881;

	double time = 0.05 * U_MYR;
	intptr_t i;
	for(i = 0; i < 10; i++) {
		double xHI = 0, xHII = 0;
		lambdaHI_to_xHI_xHII(step.lambdaHI, xHI, xHII);
		printf("%g %g %g %g %g %g \n", i * time,  step.ie, step.lambdaHI, xHI, xHII, step.dyGrec);

		step.heat = 0.0;
		step.dyGrec = 0.0;
		step.time = time;
		int code = step_evolve (&step);
		if(!code) {
			WARNING("evolve failed");
		}
	}

	return 0;
}
