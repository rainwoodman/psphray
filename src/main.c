#include <stdio.h>
#include <stdint.h>
#include <gsl/gsl_errno.h>
#include "config.h"
#include "psystem.h"

extern void init();
extern void rt_switch_epoch(int i);
extern void psystem_switch_epoch(int i);
extern void run_epoch();

int main(int argc, char* argv[]) {
	cfg_init(argv[1]);
	cfg_dump_stream(stdout);

	gsl_set_error_handler_off();
	init();
	int i;
	for(i = 0; i < N_EPOCHS; i++) {
		psystem_switch_epoch(i);

		if(psys.epoch->duration <= 0.0) continue;

		rt_switch_epoch(i);

		run_epoch();
	}
	return 0;

}
