#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <gsl/gsl_errno.h>
#include "config.h"
#include "psystem.h"

extern void init();
extern void rt_switch_epoch(int i);
extern void psystem_switch_epoch(int i);
extern void run_epoch();

int main(int argc, char* argv[]) {
	int c;
	int start_epoch = 0;
	int restart_snap = -1;
	while(-1 != (c = getopt(argc, argv, "E:S:hH"))) {
		switch(c) {
			case 'E':
				start_epoch = atoi(optarg);
			break;
			case 'S':
				restart_snap= atoi(optarg);
			break;
			case 'h':
			case 'H':
			default:
				printf("%s [-E start_epoch_id] [-S restart_snapshot_id] configfile\n", argv[0]);
				return 0;
		}
	}

	cfg_init(argv[optind]);
	cfg_dump_stream(stdout);

	gsl_set_error_handler_off();
	init();
	int i;
	for(i = start_epoch; i < N_EPOCHS; i++) {
		psystem_switch_epoch(i);

		if(psys.epoch->duration <= 0.0) continue;

		rt_switch_epoch(i);

		run_epoch();
	}
	return 0;

}
