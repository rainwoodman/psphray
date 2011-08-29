#include <messages.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "config.h"
#include <math.h>
#include <array.h>

#include <gsl/gsl_randist.h>

typedef struct {
	char * name;
	char * filename;
	int type;  /* 0 for multi, 1 for mono */
	ARRAY_DEFINE_S(freq, double)
	ARRAY_DEFINE_S(weight, double)
	gsl_ran_discrete_t * randist;
} Spec;

Spec * specs = NULL;

int N_SPECS = 0;

static double blackbody(double T, double freq) {
	return freq * freq * freq / (expm1(freq * U_RY_ENG / (C_BOLTZMANN * T * U_KELVIN)));
}

const int spec_get(const char * name) {
	int i;
	for(i = 0; i < N_SPECS; i++) {
		if(!strcmp(specs[i].name, name)) return i;
	}
	ERROR("spec %s unknown", name);
}

const char * spec_name(const int id) {
	return specs[id].name;
}

double spec_gen_freq(const int id) {
	switch(specs[id].type) {
	case 0:
		return specs[id].freq[gsl_ran_discrete(RNG, specs[id].randist)];
	case 1:
		return specs[id].freq[0];
	}
}

void spec_dump(int spec) {
	int i;
	for(i = 0; i < specs[spec].freq_length; i++) {
		printf("%lg %lg\n", specs[spec].freq[i], specs[spec].weight[i]);
	}
}
static void spec_from_file(Spec * newitem, const char * filename);
static void spec_from_config_setting(Spec * newitem, config_setting_t * setting);
void spec_init() {
	config_setting_t * list = config_lookup(CFG, "specs");
	if(list == NULL) {
		return;
	}
	N_SPECS = config_setting_length(list);
	specs = calloc(N_SPECS, sizeof(Spec));
	int i;
	for(i = 0; i < N_SPECS; i++) {
		config_setting_t * ele = config_setting_get_elem(list, i);
		specs[i].name = config_setting_name(ele);
		if(config_setting_type(ele) == CONFIG_TYPE_STRING) {
			specs[i].type = 0;
			spec_from_file(&specs[i], config_setting_get_string(ele));
		} else if(config_setting_type(ele) == CONFIG_TYPE_GROUP) {
			specs[i].type = 0;
			spec_from_config_setting(&specs[i], ele);
		} else {
			specs[i].type = 1;
			ARRAY_RESIZE(specs[i].freq, double, 1);
			specs[i].freq[0] = config_setting_get_float(ele);
		}
	}
	
	for(i = 0; i < N_SPECS; i++) {
		if(specs[i].type == 0) {
			specs[i].randist = gsl_ran_discrete_preproc(specs[i].weight_length, specs[i].weight);
		}
	}
}
static void spec_from_config_setting(Spec * newitem, config_setting_t * setting) {
	const char * type = NULL;
	if(!config_setting_lookup_string(setting, "type", &type)) {
		ERROR("can't parse spectra type in configfile");
	}
	newitem->type = 1;
	const int nbins = 1024;
	ARRAY_RESIZE(newitem->freq, double, nbins);
	ARRAY_RESIZE(newitem->weight, double, nbins);
	int i;
	double freqmin = 1;
	double freqmax = 16;
	double step = (freqmax - freqmin) / (nbins - 1);
	for(i = 0; i < nbins; i++) {
		newitem->freq[i] = freqmin + i * step;
	}
	if(!strcmp(type, "blackbody")) {
		double T;
		if(!config_setting_lookup_float(setting, "temperature", & T)
		&& !config_setting_lookup_float(setting, "T", &T)
		&& !config_setting_lookup_float(setting, "temp", &T)) {
			ERROR("specify temperature in blackbody spectrum");
		}
		for(i = 0; i < nbins; i++) {
			newitem->weight[i] = blackbody(T, newitem->freq[i]);
		}
	}

}
static void spec_from_file(Spec * newitem, const char * filename) {
	FILE * fp = fopen(filename, "r");
	if(fp == NULL) {
		ERROR("failed to open %s", filename);
	}
	
	int NR = 0;
	char * line = NULL;
	size_t len = 0;
	int stage = 0;

	double freq, weight;
	while(0 <= getline(&line, &len, fp)) {
		if(line[0] == '#') {
			NR++;
			continue;
		}
		switch(stage) {
		case 0:
			if(2 != sscanf(line, "%le %le", &freq, &weight)) {
				ERROR("%s format error at %d", filename, NR);
			} else {
				* ARRAY_APPEND(newitem->freq, double) = freq;
				* ARRAY_APPEND(newitem->weight, double) = weight;
			}
			break;
		case 1:
			break;
		}
		NR++;
		if(stage == 1) break;
	}
	free(line);
	fclose(fp);
}
