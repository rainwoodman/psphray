#include <stdio.h>
#include <stdint.h>
#include "config.h"
#include "psystem.h"

int main(int argc, char* argv[]) {
	cfg_init(argv[1]);
	cfg_dump_stream(stdout);

	init();
	int i;
	for(i = 0; i < N_EPOCHS; i++) {
		psystem_switch_epoch(i);

		rt_switch_epoch(i);

		run_epoch();
	}
	return 0;

}
