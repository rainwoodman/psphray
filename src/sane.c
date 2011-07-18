#include <stdio.h>
#include <stdint.h>
#include <messages.h>
#include "reader.h"
#include "config.h"
#include "psystem.h"
#include "atomic-rate.h"

extern size_t rt_trace(const float s[3], const float dir[3], const float dist, intptr_t ** pars, size_t * size);
extern int pluecker_(float dir[3], float * dist, float s2b[3], float s2t[3]);

extern PSystem psys;

int main(int argc, char* argv[]) {
	float pos[3] = {1000., 1000., 1000.};
	float s2b[3] = {-1000., -1000., -1000.};
	float s2t[3] = {1000., 1000., 1000.};
	float dir[3] = {0., 1., 0.};
	float dist = 100000.0;
	if(!pluecker_(dir, &dist, s2b, s2t)) {
		ERROR("shall be true");
	}

	ar_init("../test/atomic_rates_Hui.txt");
	MESSAGE("1 ryd verner cross section is %e", ar_verner(1.0));
	MESSAGE("1e5K HI_CI %e HII_RC_A %e HII_RCC_A %e", 
			ar_get(AR_HI_CI, 5),
			ar_get(AR_HII_RC_A, 5),
			ar_get(AR_HII_RCC_A, 5)
		);

	cfg_init("../test/configfile");
	cfg_dump("used-config");

	MESSAGE("--------unit sanity check------");
	MESSAGE("kpc /h = %f", units_parse("kpc / h"));
	MESSAGE("10^10 msun /h = %f", 1e10 * units_parse("msun / h"));
	MESSAGE("980 myr /h = %f", 980 * units_parse("myr / h"));
	MESSAGE("--------unit sanity check------");

	Reader * r = reader_new("massiveblack");
	reader_open(r, "../test/snapdir_680/snapshot_680.0");
	ReaderConstants * c = reader_constants(r);
	MESSAGE("-------reader sanity check--------");
	MESSAGE("boxsize = %f", c->boxsize);
	MESSAGE("redshift = %f", c->redshift);
	MESSAGE("Ntot[0] = %lu", c->Ntot[0]);
	MESSAGE("N[0] = %lu", c->N[0]);
	MESSAGE("Nfiles = %lu", c->Nfiles);
	MESSAGE("h = %f", c->h);
	MESSAGE("-------reader sanity check--------");

	float *ie = reader_alloc(r, "ie", 0);
	float *ye = reader_alloc(r, "ye", 0);
	reader_read(r, "ie", 0, ie);
	reader_read(r, "ye", 0, ye);

	float Tmin = 1000000000;
	float Tmax = -1;
	intptr_t i;
	for(i = 0; i < reader_npar(r, 0); i++) {
		float T = ieye2T(ie[i], ye[i]);
		if(T > Tmax) Tmax = T;
		if(T < Tmin) Tmin = T;
	}
	free(ie);
	free(ye);
	reader_destroy(r);
	printf("%f %f\n", Tmax, Tmin);

	init();
	psystem_switch_epoch(0);

	rt_switch_epoch(0);

	run();

#if 0
	intptr_t * ipars = NULL;
	size_t size =0;
	size_t length = rt_trace(pos, dir, dist, &ipars, &size);
	for(i = 0; i < length; i++) {
		intptr_t ipar = ipars[i];
		printf("%g %g %g\n", psys.pos[ipar][0], psys.pos[ipar][1], psys.pos[ipar][2]);
	}
	free(ipars);
	return 0;
	psystem_switch_epoch(1);
	rt_switch_epoch(1);
#endif

	return 0;
}
