#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <messages.h>
#include <libconfig/libconfig.h>
#include <gsl/gsl_rng.h>

extern void units_init();
extern void spec_init();
extern void epochs_init();
extern double units_parse(const char * units);
extern void ar_init(const char * filename);

config_t CFG[] = {{0}};

gsl_rng * RNG = NULL;

int CFG_WRITE_INIT = 0;
int CFG_ISOTHERMAL = 0;
int CFG_ADIABATIC = 0;
int CFG_DISABLE_2ND_GEN_PHOTONS = 0;
int CFG_COMOVING = 1;

config_setting_t * config_ensure(config_t * config, char * path, int type);
#define config_ensure_int(c, p, v)  if(!config_lookup(c, p)) config_setting_set_int(config_ensure(c, p, CONFIG_TYPE_INT), v)
#define config_ensure_int64(c, p, v)  if(!config_lookup(c,p)) config_setting_set_int64(config_ensure(c, p, CONFIG_TYPE_INT64), v)
#define config_ensure_float(c, p, v)  if(!config_lookup(c,p)) config_setting_set_float(config_ensure(c, p, CONFIG_TYPE_FLOAT), v)
#define config_ensure_bool(c, p, v)  if(!config_lookup(c,p)) config_setting_set_bool(config_ensure(c, p, CONFIG_TYPE_BOOL), v)
#define config_ensure_string(c, p, v)  if(!config_lookup(c,p)) config_setting_set_string(config_ensure(c, p, CONFIG_TYPE_STRING), v)

void cfg_init(char * filename) {
	config_init(CFG);
	config_set_auto_convert(CFG, CONFIG_TRUE);
	if(!config_read_file(CFG, filename)) {
		ERROR("%s: %d: %s", config_error_file(CFG), config_error_line(CFG), config_error_text(CFG));
	}
	config_ensure        (CFG, "psphray", CONFIG_TYPE_GROUP);
	config_ensure_string (CFG, "psphray.atomicRates", "<filename>");
	config_ensure_int64  (CFG, "psphray.seed", 123456);
	config_ensure_bool   (CFG, "psphray.writeInit", CONFIG_TRUE);
	config_ensure_bool   (CFG, "psphray.disable2ndGenPhotons", CONFIG_FALSE);
	config_ensure_bool   (CFG, "psphray.isothermal", CONFIG_TRUE);
	config_ensure_bool   (CFG, "psphray.adiabatic", CONFIG_FALSE);

	config_ensure        (CFG, "cosmology", CONFIG_TYPE_GROUP);
	config_ensure_bool   (CFG, "cosmology.comoving", CONFIG_TRUE);
	config_ensure_float  (CFG, "cosmology.omegaL", 0.74);
	config_ensure_float  (CFG, "cosmology.omegaM", 0.26);
	config_ensure_float  (CFG, "cosmology.omegaB", 0.044);
	config_ensure_float  (CFG, "cosmology.h", 0.72);
	config_ensure_float  (CFG, "cosmology.hmf", 0.76);
	config_ensure        (CFG, "box", CONFIG_TYPE_GROUP);
	config_ensure_float  (CFG, "box.boxsize", -1.0);
	config_ensure_string (CFG, "box.boundary", "vaccum");
	config_ensure        (CFG, "units", CONFIG_TYPE_GROUP);
	config_ensure_float  (CFG, "units.massGramh", 1.9847005219450705e+43);
	config_ensure_float  (CFG, "units.lengthCMh", 3.0835455558318480e+21);
	config_ensure_float  (CFG, "units.timeSh", 3.08568025e+16);

	config_lookup_bool(CFG, "psphray.writeInit", &CFG_WRITE_INIT);
	config_lookup_bool(CFG, "psphray.disable2ndGenPhotons", &CFG_DISABLE_2ND_GEN_PHOTONS);
	config_lookup_bool(CFG, "psphray.isothermal", &CFG_ISOTHERMAL);
	config_lookup_bool(CFG, "psphray.adiabatic", &CFG_ADIABATIC);
	config_lookup_bool(CFG, "cosmology.comoving", &CFG_COMOVING);

	int64_t seed = 123456;
	config_lookup_int64(CFG, "psphray.seed", &seed);
	RNG = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(RNG, seed);

	units_init();
	const char * arfilename = NULL;
	config_lookup_string(CFG, "psphray.atomicRates", &arfilename);
	ar_init(arfilename);

	spec_init();

	epochs_init();
}

void cfg_dump(char * filename) {
	config_write_file(CFG, filename);
}

void cfg_dump_stream(FILE * file) {
	config_write(CFG, file);
}

config_setting_t * config_ensure(config_t * config, char * path, int type) {
	config_setting_t * st = config_lookup(config, path);
	if(st == NULL) {
		char * pntpath = strdup(path);
		char * dot = rindex(pntpath, '.');
		config_setting_t * pnt = config_root_setting(config);
		if(dot != NULL) {
			*dot = 0;
			char * name = dot + 1;
			pnt = config_lookup(CFG, pntpath);
			if(pnt == NULL) {
				ERROR("RUNTIME: %s not found", pntpath);
			}
			st = config_setting_add(pnt, name, type);
		}
		free(pntpath);
	}
	return st;
}

config_setting_t * config_setting_ensure_member(config_setting_t * e, char * member, int type) {
	config_setting_t * mem = config_setting_get_member(e, member);
	if(mem == NULL) {
		mem = config_setting_add(e, member, type);
	}
	return mem;
}

double config_setting_parse_units(config_setting_t * e) {
	if(config_setting_is_list(e)) {
		double value = config_setting_get_float_elem(e, 0);
		const char * units = config_setting_get_string_elem(e, 1);
		return value * units_parse(units);
	} else {
		double value = config_setting_get_float(e);
		return value;
	}
}
