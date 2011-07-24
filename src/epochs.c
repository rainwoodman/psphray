#include <messages.h>
#include <stdlib.h>
#include <stdint.h>
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

static void epochs_input(config_setting_t * list, int i) {
	config_setting_t * e = config_setting_get_elem(list, i);
	config_setting_t * s = config_setting_get_member(e, "snapshot");
	const char * snapshot = config_setting_get_string_elem(s, 0);
	const char * format = config_setting_get_string_elem(s, 1);
	const char * source = NULL;
	double redshift;
	size_t ngas = 0;
	int nticks = 0;
	int nray = 1;

	config_setting_lookup_string(e, "source", &source);
	config_setting_lookup_int(e, "nticks", &nticks);
	config_setting_lookup_int(e, "nray", &nray);


	Reader * r = reader_new(format);
	char * fname = reader_make_filename(snapshot, 0);
	reader_open(r, fname);
	free(fname);
	ReaderConstants * c = reader_constants(r);

	redshift = c->redshift;
	ngas = c->Ntot[0];

	reader_destroy(r);

	EPOCHS[i].snapshot = snapshot;
	EPOCHS[i].format = format;
	EPOCHS[i].source = source;
	EPOCHS[i].redshift = redshift;
	EPOCHS[i].age = z2t(redshift);
	EPOCHS[i].ngas = ngas;
	EPOCHS[i].nray = nray;
	EPOCHS[i].nticks = nticks;
}

static void epochs_duration(config_setting_t * list, int i) {
	config_setting_t * e = config_setting_get_elem(list, i);
	config_setting_t * d = config_setting_get_member(e, "duration");
	if(d) {
		EPOCHS[i].duration = config_setting_parse_units(d);
		return;
	} 
	if(i != N_EPOCHS - 1) {
		EPOCHS[i].duration = z2t(EPOCHS[i+1].redshift) - z2t(EPOCHS[i].redshift);
		return;
	}
	EPOCHS[i].duration = 0;
}

static void epochs_output(config_setting_t * list, int i) {
	config_setting_t * e = config_setting_get_elem(list, i);
	const char * output = NULL;
	int output_nfiles = 0;
	int output_nsteps = 0;
	intptr_t * output_steps = NULL;
	double duration = EPOCHS[i].duration;
	size_t nticks =EPOCHS[i].nticks;

	config_setting_t * o = config_setting_get_member(e, "output");
	if(o == NULL) {
		MESSAGE("EPOCH %d(%s) has no output", i, EPOCHS[i].snapshot);
		return;
	}

	config_setting_lookup_string(o, "filename", &output);
	config_setting_lookup_int(o, "nfiles", &output_nfiles);

	config_setting_t * s = config_setting_get_member(o, "steps");
	
	int j;
	switch(config_setting_type(s)){
	case CONFIG_TYPE_INT:
	case CONFIG_TYPE_INT64:
	case CONFIG_TYPE_FLOAT:
		output_nsteps = config_setting_get_int(s);
		output_steps = calloc(sizeof(size_t), output_nsteps);
		for(j = 0; j < output_nsteps; j++) {
			output_steps[j] = nticks * (j+1)/ output_nsteps;
		}
	break;
	case CONFIG_TYPE_GROUP:
	break;
	case CONFIG_TYPE_LIST:
		output_nsteps = config_setting_length(s);
		output_steps = calloc(sizeof(size_t), output_nsteps);
		for(j = 0; j < output_nsteps; j++) {
			double time = 0.0;
			config_setting_t * st = config_setting_get_elem(s, j);
			time = config_setting_parse_units(st);
			output_steps[j] = time * nticks / duration + 0.5; /*rounding*/
			if(output_steps[j] > nticks) output_steps[j] = nticks;
		}
	break;
	}

	EPOCHS[i].output.filename = output;
	EPOCHS[i].output.nfiles = output_nfiles;
	EPOCHS[i].output.steps = output_steps;
	EPOCHS[i].output.nsteps = output_nsteps;

}

void epochs_init() {
	config_setting_t * list = config_lookup(CFG, "epochs");
	if(list == NULL) {
		return;
	}
	N_EPOCHS = config_setting_length(list);
	int i;
	EPOCHS = calloc(sizeof(Epoch), N_EPOCHS);
	for(i = 0; i < N_EPOCHS; i++) {
		epochs_input(list, i);
	}
	for(i = 0; i < N_EPOCHS; i++) {
		epochs_duration(list, i);
	}
	for(i = 0; i < N_EPOCHS; i++) {
		epochs_output(list, i);
	}

	for(i = 0; i < N_EPOCHS; i++) {
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
		intptr_t j;
		for(j = 0; j < EPOCHS[i].output.nsteps; j++)
			MESSAGE("epoch %d, output %d @ [tick=%ld, time=%g myr]",
				i, j + 1, EPOCHS[i].output.steps[j],
				units_format(EPOCHS[i].output.steps[j] * EPOCHS[i].duration / EPOCHS[i].nticks, "myr")
				);
	}
}

