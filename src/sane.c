#include <stdio.h>
#include <stdint.h>
#include <messages.h>
#include "reader.h"
#include "config.h"
#include "psystem.h"

extern void spec_dump(int spec);
extern void solve_u_v(double * d, double * u, double * v);
extern double sph_depth(const double r_h);
int main(int argc, char* argv[]) {
	intptr_t i;
	cfg_init(argv[1]);
	cfg_dump_stream(stdout);

	MESSAGE("HMF = %g", C_HMF);
	MESSAGE("1 ryd verner cross section is %e", xs_get(XS_HI, 1.0));
	MESSAGE("1.81 ryd HeI cross section is %e", xs_get(XS_HEI, 1.81));
	MESSAGE("1 ryd energy = %e", U_RY_ENG);
	MESSAGE("1 Kelvin energy = %e", C_BOLTZMANN * 1 * U_KELVIN);
	MESSAGE("1e4K HI_CI %e HII_RC_A %e HII_RCC_A %e", 
			ar_get(AR_HI_CI, 4),
			ar_get(AR_HII_RC_A, 4),
			ar_get(AR_HII_RCC_A, 4)
		);
	MESSAGE("1e4K HeI_CI %e HeII_RC_A %e HeII_RC_B %e", 
			ar_get(AR_HEI_CI, 4),
			ar_get(AR_HEII_RC_A, 4),
			ar_get(AR_HEII_RC_B, 4)
		);
	MESSAGE("1e4K HeII_CI %e HeIII_RC_A %e HeIII_RC_B %e", 
			ar_get(AR_HEII_CI, 4),
			ar_get(AR_HEIII_RC_A, 4),
			ar_get(AR_HEIII_RC_B, 4)
		);

	MESSAGE("sph depth %g", sph_depth(0.66));
	MESSAGE("10 < 1 = %d", (1. > 10.) - (10. > 1.));
	MESSAGE("--------unit sanity check------");
	MESSAGE("kpc /h = %f", units_parse("kpc / h"));
	MESSAGE("10^10 msun /h = %f", 1e10 * units_parse("msun / h"));
	MESSAGE("980 myr /h = %f", 980 * units_parse("myr / h"));
	MESSAGE("density 1 = %f proton / cm^3", 1 / units_parse(" mproton / cm cm cm"));
	MESSAGE("--------unit sanity check------");

	for(i = 0; i < 1000; i++) {
		printf("%g %g\n", i * 0.01, xs_get(XS_HEI, i * 0.01));
	}
#if 0
	double d[] = {1.0, 0.0, 0.0};
	double u[3];
	double v[3];
	solve_u_v(d, u, v);
	MESSAGE("d = (%g %g %g), u=(%g %g %g), v=(%g %g %g)",
		d[0],d[1],d[2],
		u[0],u[1],u[2],
		v[0],v[1],v[2]);

	const int spec = spec_get("sun");
	spec_dump(spec);
	for(i = 0; i < 100; i++) {
		printf("%lg\n", spec_gen_freq(spec));
	}
#endif
	return 0;
}
