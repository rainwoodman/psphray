#include <stdio.h>
#include <stdint.h>
#include <messages.h>
#include "reader.h"
#include "config.h"

int main(int argc, char* argv[]) {
	cfg_init("../tests/configfile");
	cfg_dump("used-config");

	MESSAGE("--------unit sanity check------");
	MESSAGE("kpc /h = %f", units_parse("kpc / h"));
	MESSAGE("10^10 msun /h = %f", 1e10 * units_parse("msun / h"));
	MESSAGE("980 myr /h = %f", 980 * units_parse("myr / h"));
	MESSAGE("--------unit sanity check------");

	Reader * r = reader_new("massiveblack");
	reader_open(r, "../test/snapshot_036.0");
	float (*pos)[3] = reader_alloc(r, "pos", -1);
	size_t l = reader_length(r, "pos");

	reader_read(r, "pos", -1, pos);
	ReaderConstants * c = reader_constants(r);
	MESSAGE("-------reader sanity check--------");
	MESSAGE("boxsize = %f", c->boxsize);
	MESSAGE("redshift = %f", c->redshift);
	MESSAGE("Ntot[0] = %lu", c->Ntot[0]);
	MESSAGE("N[0] = %lu", c->N[0]);
	MESSAGE("Nfiles = %lu", c->Nfiles);
	MESSAGE("h = %f", c->h);
	MESSAGE("-------reader sanity check--------");

	intptr_t i;
	float x, y, z;
	for(i = 0; i < l; i++) {
		x += pos[i][0];
		y += pos[i][1];
		z += pos[i][2];
	}
	printf("%f %f %f\n", x/l, y/l, z/l);
	return 0;
}
