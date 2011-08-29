#include <stdio.h>
#include <stdint.h>
#include <messages.h>
#include "reader.h"
#include "config.h"
#include "psystem.h"

extern void spec_dump(int spec);
int main(int argc, char* argv[]) {
	intptr_t i;
	cfg_init(argv[1]);
//	cfg_dump_stream(stdout);

	MESSAGE("HMF = %g", C_HMF);
	MESSAGE("1 ryd verner cross section is %e", xs_get(XS_HI, 1.0));
	MESSAGE("1 ryd energy = %e", U_RY_ENG);
	MESSAGE("1 Kelvin energy = %e", C_BOLTZMANN * 1 * U_KELVIN);
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

	const int spec = spec_get("sun");
	spec_dump(spec);
	return 0;
}
