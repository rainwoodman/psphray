#include <messages.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "reader.h"

Epoch * EPOCHS = NULL;
int N_EPOCHS = 0;

#define config_setting_ensure_int(e, m, v)  if(!config_setting_get_member(e, m)) config_setting_set_int(config_setting_ensure_member(c, m, CONFIG_TYPE_INT), v)
#define config_setting_ensure_int64(e, m, v)  if(!config_setting_get_member(e, m)) config_setting_set_int64(config_setting_ensure_member(c, m, CONFIG_TYPE_INT64), v)
#define config_setting_ensure_float(e, m, v)  if(!config_setting_get_member(e, m)) config_setting_set_float(config_setting_ensure_member(c, m, CONFIG_TYPE_FLOAT), v)
#define config_setting_ensure_string(e, m, v)  if(!config_setting_get_member(e, m)) config_setting_set_string(config_setting_ensure_member(c, m, CONFIG_TYPE_STRING), v)
#define config_setting_ensure_bool(e, m, v)  if(!config_setting_get_member(e, m)) config_setting_set_bool(config_setting_ensure_member(c, m, CONFIG_TYPE_BOOL), v)

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
		const char * output = NULL;
		int output_nfiles = 0;
		int output_nsteps = 0;
		double duration = -1;
		double redshift;
		size_t ngas = 0;
		int nticks = 0;
		int nray = 1;

		config_setting_lookup_string(e, "source", &source);
		config_setting_lookup_int(e, "nticks", &nticks);
		config_setting_lookup_int(e, "nray", &nray);

		config_setting_t * o = config_setting_get_member(e, "output");
		if(o != NULL) {
			config_setting_lookup_string(o, "filename", &output);
			config_setting_lookup_int(o, "nsteps", &output_nsteps);
			config_setting_lookup_int(o, "nfiles", &output_nfiles);
		} else {
			MESSAGE("EPOCH %d(%s) has no output", i, snapshot);
		}

		Reader * r = reader_new(format);
		char * fname = reader_make_filename(snapshot, 0);
		reader_open(r, fname);
		free(fname);
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
		EPOCHS[i].nray = nray;
		EPOCHS[i].nticks = nticks;
		EPOCHS[i].output.filename = output;
		EPOCHS[i].output.nfiles = output_nfiles;
		EPOCHS[i].output.nsteps = output_nsteps;
	}
	for(i = 0; i < N_EPOCHS; i++) {
		if(EPOCHS[i].duration == -1 && i != N_EPOCHS - 1) {
			EPOCHS[i].duration = z2t(EPOCHS[i+1].redshift) - z2t(EPOCHS[i].redshift);
		}
		MESSAGE("epoch %d, redshift %g, snapshot %s, format %s, source %s, ticks %u, ngas %lu, age %f [myr], duration %f [myr]",
		i, 
		EPOCHS[i].redshift,
		EPOCHS[i].snapshot,
		EPOCHS[i].format,
		EPOCHS[i].source,
		EPOCHS[i].nticks,
		EPOCHS[i].ngas,
		units_format(EPOCHS[i].age, "myr"),
		units_format(EPOCHS[i].duration, "myr"));
		int j;
		EPOCHS[i].output.steps = calloc(sizeof(size_t), EPOCHS[i].output.nsteps);
		for(j = 0; j < EPOCHS[i].output.nsteps; j++) {
			EPOCHS[i].output.steps[j] = EPOCHS[i].nticks * (j+1)/ EPOCHS[i].output.nsteps;
		}
	}
}

