#include <stdio.h>
#include <stdint.h>
#include "config.h"
#include "psystem.h"

extern void psystem_switch_epoch(int i);

int main(int argc, char* argv[]) {
	cfg_init(argv[1]);
	cfg_dump_stream(stdout);

	psystem_switch_epoch(0);

	intptr_t ipar;
	intptr_t start = 0;
	intptr_t end = psys.npar / 100;
	for(ipar = start; ipar < end; ipar++) {
		printf("%g %g %g\n", psys.pos[ipar][0], psys.pos[ipar][1], psys.pos[ipar][2]);
	}
	return 0;

}
