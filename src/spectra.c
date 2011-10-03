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
	ARRAY_DEFINE_S(eng, double)
	ARRAY_DEFINE_S(lum, double)
	ARRAY_DEFINE_S(lumwt, double)
	ARRAY_DEFINE_S(Nwt, double)
	gsl_ran_discrete_t * randist;
	double band_min;
	double band_max;
	double N_lum;
} Spec;

Spec * specs = NULL;

int N_SPECS = 0;

static double blackbody(double T, double eng) {
	return eng * eng * eng / (expm1(eng / (C_BOLTZMANN * T * U_KELVIN)));
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

/* returns the energy of the photon, not the freq! */
double spec_gen_freq(const int id) {
	switch(specs[id].type) {
	case 0:
		return specs[id].eng[gsl_ran_discrete(RNG, specs[id].randist)];
	case 1:
		return specs[id].eng[0];
	}
}

double spec_N_from_lum(const int id, double lum) {
	return lum * specs[id].N_lum;
}
double spec_lum_from_N(const int id, double N) {
	return N / specs[id].N_lum;
}

void spec_dump(int spec) {
	int i;
	for(i = 0; i < specs[spec].eng_length; i++) {
		printf("%lg %lg\n", specs[spec].eng[i], specs[spec].lum[i]);
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
			ARRAY_RESIZE(specs[i].eng, double, 1);
			ARRAY_RESIZE(specs[i].lum, double, 1);
			specs[i].eng[0] = config_setting_get_float(ele);
			specs[i].lum[0] = 1.0;
		}
	}
	
	for(i = 0; i < N_SPECS; i++) {
		if(specs[i].type != 0) {
			specs[i].N_lum = 1.0 / specs[i].eng[0];
			continue;
		}
		double lumwtsum = 0.0;
		double Nwtsum = 0.0;
		int j;
		ARRAY_RESIZE(specs[i].Nwt, double, specs[i].lum_length);
		ARRAY_RESIZE(specs[i].lumwt, double, specs[i].lum_length);
		for(j = 0; j < specs[i].Nwt_length; j++) {
			specs[i].lumwt[j] = specs[i].lum[j];
			specs[i].Nwt[j] = specs[i].lum[j] / specs[i].eng[j];
			if(j == 0) {
				specs[i].lumwt[j] *= 2 * (specs[i].eng[j+1] - specs[i].eng[j]);
				specs[i].Nwt[j] *= 2 * (specs[i].eng[j+1] - specs[i].eng[j]);
			} else if(j == specs[i].lumwt_length - 1) {
				specs[i].lumwt[j] *= 2 * (specs[i].eng[j] - specs[i].eng[j-1]);
				specs[i].Nwt[j] *= 2 * (specs[i].eng[j] - specs[i].eng[j-1]);
			} else {
				specs[i].lumwt[j] *= (specs[i].eng[j + 1] - specs[i].eng[j-1]);
				specs[i].Nwt[j] *= (specs[i].eng[j + 1] - specs[i].eng[j-1]);
			}
			lumwtsum += specs[i].lumwt[j];
			Nwtsum += specs[i].Nwt[j];
		}
		specs[i].N_lum = Nwtsum / lumwtsum;
		specs[i].band_min = specs[i].eng[0];
		specs[i].band_max  = specs[i].eng[specs[i].eng_length - 1];
		specs[i].randist = gsl_ran_discrete_preproc(specs[i].Nwt_length, specs[i].Nwt);
	}
}

static void spec_from_config_setting(Spec * newitem, config_setting_t * setting) {
	const char * type = NULL;
	if(!config_setting_lookup_string(setting, "type", &type)) {
		ERROR("can't parse spectra type in configfile");
	}
	double engmin = 0;
	double engmax = 0;
	if(!config_setting_parse_units_member(setting, "min", &engmin)) {
		ERROR("can't parse engmin in configfile");
	}
	if(!config_setting_parse_units_member(setting, "max", &engmax)) {
		ERROR("can't parse max in configfile");
	}
	const int nbins = 1024;
	ARRAY_RESIZE(newitem->eng, double, nbins);
	ARRAY_RESIZE(newitem->lum, double, nbins);
	newitem->band_min = engmin;
	newitem->band_max = engmax;
	int i;
	double step = (engmax - engmin) / (nbins - 1);
	for(i = 0; i < nbins; i++) {
		newitem->eng[i] = engmin + i * step;
	}
	if(!strcmp(type, "blackbody")) {
		double T;
		if(!config_setting_lookup_float(setting, "temperature", & T)
		&& !config_setting_lookup_float(setting, "T", &T)
		&& !config_setting_lookup_float(setting, "temp", &T)) {
			ERROR("specify temperature in blackbody spectrum");
		}
		for(i = 0; i < nbins; i++) {
			newitem->lum[i] = blackbody(T, newitem->eng[i]);
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

	double freq, N;
	while(0 <= getline(&line, &len, fp)) {
		if(line[0] == '#') {
			NR++;
			continue;
		}
		switch(stage) {
		case 0:
			if(2 != sscanf(line, "%le %le", &freq, &N)) {
				ERROR("%s format error at %d", filename, NR);
			} else {
				* ARRAY_APPEND(newitem->eng, double) = freq * U_RY_ENG;
				* ARRAY_APPEND(newitem->lum, double) = N * freq * U_RY_ENG;
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
