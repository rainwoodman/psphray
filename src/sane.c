#include <stdio.h>
#include <stdint.h>
#include <messages.h>
#include "reader.h"
#include "config.h"
#include "psystem.h"

extern size_t rt_trace(const float s[3], const float dir[3], const float dist, intptr_t ** pars, size_t * size);
extern int pluecker_(float dir[3], float * dist, float s2b[3], float s2t[3]);

extern PSystem psys;

int main(int argc, char* argv[]) {
	intptr_t i;
	cfg_init("../test/uniform.config");
	cfg_dump("used-config");

	MESSAGE("1 ryd verner cross section is %e", ar_verner(1.0));
	MESSAGE("1e5K HI_CI %e HII_RC_A %e HII_RCC_A %e", 
			ar_get(AR_HI_CI, 5),
			ar_get(AR_HII_RC_A, 5),
			ar_get(AR_HII_RCC_A, 5)
		);

	MESSAGE("--------unit sanity check------");
	MESSAGE("kpc /h = %f", units_parse("kpc / h"));
	MESSAGE("10^10 msun /h = %f", 1e10 * units_parse("msun / h"));
	MESSAGE("980 myr /h = %f", 980 * units_parse("myr / h"));
	MESSAGE("density 1 = %f proton / cm^3", 1 / units_parse(" mproton / cm cm cm"));
	MESSAGE("--------unit sanity check------");

	init();
	psystem_switch_epoch(0);

	rt_switch_epoch(0);

	run();
	return 0;

	psystem_switch_epoch(1);
	rt_switch_epoch(1);

	return 0;
}
