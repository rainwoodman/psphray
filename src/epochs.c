#include <messages.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "reader.h"

Epoch * EPOCHS = NULL;
int N_EPOCHS = 0;
void epochs_init() {
	config_setting_t * list = config_lookup(CFG, "epochs");
	N_EPOCHS = config_setting_length(list);
	int i;
	EPOCHS = calloc(sizeof(Epoch), N_EPOCHS);
	for(i = 0; i < N_EPOCHS; i++) {
		config_setting_t * e = config_setting_get_elem(list, i);
		config_setting_t * s = config_setting_get_member(e, "snapshot");
		const char * snapshot = config_setting_get_string_elem(s, 0);
		const char * format = config_setting_get_string_elem(s, 1);
		const char * source = NULL;
		double duration = -1;
		double redshift;
		size_t ngas = 0;

		int nrays = 0;
		config_setting_lookup_string(e, "source", &source);
		config_setting_lookup_int(e, "nrays", &nrays);

		Reader * r = reader_new(format);
		reader_open(r, snapshot);
		ReaderConstants * c = reader_constants(r);

		redshift = c->redshift;
		ngas = c->Ntot[0];

		reader_destroy(r);

		config_setting_t * d = config_setting_get_member(e, "duration");
		if(d) {
			duration = config_setting_parse_units(d);
		}
		EPOCHS[i].snapshot = snapshot;
		EPOCHS[i].format = format;
		EPOCHS[i].source = source;
		EPOCHS[i].redshift = redshift;
		EPOCHS[i].age = z2t(redshift);
		EPOCHS[i].duration = duration;
		EPOCHS[i].ngas = ngas;
		EPOCHS[i].nrays = nrays;

	}
	for(i = 0; i < N_EPOCHS; i++) {
		if(EPOCHS[i].duration == -1 && i != N_EPOCHS - 1) {
			EPOCHS[i].duration = z2t(EPOCHS[i+1].redshift) - z2t(EPOCHS[i].redshift);
		}
		MESSAGE("epoch %d, redshift %g, snapshot %s, format %s, source %s, nrays %u, ngas %lu, age %f [myr], duration %f [myr]",
		i, 
		EPOCHS[i].redshift,
		EPOCHS[i].snapshot,
		EPOCHS[i].format,
		EPOCHS[i].source,
		EPOCHS[i].nrays,
		EPOCHS[i].ngas,
		units_format(EPOCHS[i].age, "myr"),
		units_format(EPOCHS[i].duration, "myr"));
	}
}
