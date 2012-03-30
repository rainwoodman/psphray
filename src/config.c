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
extern void xs_init(const char * filename);
extern void eos_init(const char * filename);
extern void lte_init();
extern void secion_init();

config_t CFG[] = {{0}};

gsl_rng * RNG = NULL;

int CFG_WRITE_INIT = CONFIG_TRUE;
int CFG_ISOTHERMAL = CONFIG_TRUE;
int CFG_ADIABATIC = CONFIG_FALSE;
int CFG_ON_THE_SPOT = CONFIG_FALSE;
int CFG_COMOVING = CONFIG_TRUE;
int CFG_DUMP_HOTSPOTS = CONFIG_FALSE;
int CFG_H_ONLY = CONFIG_FALSE;
int CFG_TRACE_ONLY = CONFIG_FALSE;
int CFG_NO_PHOTON = CONFIG_FALSE;
int CFG_SECONDARY_IONIZATION = CONFIG_TRUE;
int CFG_OCTTREE_PPC = 8;
int CFG_DISABLE_LTE = CONFIG_FALSE;
int CFG_PRINT_RAYS = CONFIG_FALSE;
int CFG_ENABLE_EOS = CONFIG_FALSE;
double CFG_EOS_CLOUD_TEMPERATURE = 1e3;
double CFG_RECOMBINE_THRESHOLD = 0.1;
double CFG_BOOST_SOURCE_FACTOR = 1.0;

double CFG_FAKE_TEMPERATURE = -1;

int64_t CFG_SEED = 123456;
double C_1_CMH = 3.0835455558318480e+21;
double C_1_GRAMH = 1.9847005219450705e+43;
double C_1_SECH = 3.08568025e+16;
double C_OMEGA_L = 0.74;
double C_OMEGA_M = 0.26;
double C_OMEGA_B = 0.044;
double C_H = 0.72;
double C_HMF = 0.76;
double C_HEMF = 0.24;

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
	config_ensure_string (CFG, "psphray.atomicRates", NULL);
	config_ensure_string (CFG, "psphray.crossSections", NULL);
	config_ensure_int64  (CFG, "psphray.seed", CFG_SEED);
	config_ensure_bool   (CFG, "psphray.writeInit", CFG_WRITE_INIT);
	config_ensure_bool   (CFG, "psphray.onTheSpot", CFG_ON_THE_SPOT);
	config_ensure_bool   (CFG, "psphray.isothermal", CFG_ISOTHERMAL);
	config_ensure_bool   (CFG, "psphray.HOnly", CFG_H_ONLY);
	config_ensure_bool   (CFG, "psphray.adiabatic", CFG_ADIABATIC);
	config_ensure_bool   (CFG, "psphray.dumpHotspots", CFG_DUMP_HOTSPOTS);
	config_ensure_bool   (CFG, "psphray.traceOnly", CFG_TRACE_ONLY);
	config_ensure_bool   (CFG, "psphray.noPhoton", CFG_NO_PHOTON);
	config_ensure_bool   (CFG, "psphray.disableLTE", CFG_DISABLE_LTE);
	config_ensure_bool   (CFG, "psphray.printRays", CFG_PRINT_RAYS);
	config_ensure_bool   (CFG, "psphray.secondaryIonization", CFG_SECONDARY_IONIZATION);
	config_ensure        (CFG, "psphray.octtree", CONFIG_TYPE_GROUP);
	config_ensure_int64  (CFG, "psphray.octtree.particlesPerCell", CFG_OCTTREE_PPC);
	config_ensure_float  (CFG, "psphray.fakeTemperature", CFG_FAKE_TEMPERATURE);
	config_ensure_float  (CFG, "psphray.recombineThreshold", CFG_RECOMBINE_THRESHOLD);
	config_ensure_float  (CFG, "psphray.boostSourceFactor", CFG_BOOST_SOURCE_FACTOR);

	config_ensure        (CFG, "psphray.eos", CONFIG_TYPE_GROUP);
	config_ensure_string (CFG, "psphray.eos.filename", NULL);
	config_ensure_float (CFG, "psphray.eos.cloudTemperature", CFG_EOS_CLOUD_TEMPERATURE);

	config_ensure        (CFG, "cosmology", CONFIG_TYPE_GROUP);
	config_ensure_bool   (CFG, "cosmology.comoving", CFG_COMOVING);
	config_ensure_float  (CFG, "cosmology.omegaL", C_OMEGA_L);
	config_ensure_float  (CFG, "cosmology.omegaM", C_OMEGA_M);
	config_ensure_float  (CFG, "cosmology.omegaB", C_OMEGA_B);
	config_ensure_float  (CFG, "cosmology.h", C_H);
	config_ensure_float  (CFG, "cosmology.hmf", C_HMF);
	config_ensure        (CFG, "box", CONFIG_TYPE_GROUP);
	config_ensure_float  (CFG, "box.boxsize", -1.0);
	config_ensure_string (CFG, "box.boundary", "vaccum");
	config_ensure        (CFG, "units", CONFIG_TYPE_GROUP);
	config_ensure_float  (CFG, "units.massGramh", C_1_GRAMH);
	config_ensure_float  (CFG, "units.lengthCMh", C_1_CMH);
	config_ensure_float  (CFG, "units.timeSh", C_1_SECH);

	config_lookup_int64(CFG, "psphray.seed", &CFG_SEED);
	RNG = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(RNG, CFG_SEED);

	/* init the units first */ 
	config_lookup_float(CFG, "cosmology.h", &C_H);
	config_lookup_float(CFG, "cosmology.hmf", &C_HMF);
	C_HEMF = 1 - C_HMF;
	config_ensure_float(CFG, "cosmology.hemf", C_HEMF);
	config_lookup_float(CFG, "cosmology.hemf", &C_HEMF);
	config_lookup_float(CFG, "cosmology.omegaL", &C_OMEGA_L);
	config_lookup_float(CFG, "cosmology.omegaM", &C_OMEGA_M);
	config_lookup_float(CFG, "cosmology.omegaB", &C_OMEGA_B);
	
	config_lookup_float(CFG, "units.lengthCMh", &C_1_CMH);
	config_lookup_float(CFG, "units.massGramh", &C_1_GRAMH);
	config_lookup_float(CFG, "units.timeSh", &C_1_SECH);

	units_init();

	/* then config the run */
	config_lookup_bool(CFG, "psphray.writeInit", &CFG_WRITE_INIT);
	config_lookup_bool(CFG, "psphray.onTheSpot", &CFG_ON_THE_SPOT);
	config_lookup_bool(CFG, "psphray.isothermal", &CFG_ISOTHERMAL);
	config_lookup_bool(CFG, "psphray.adiabatic", &CFG_ADIABATIC);
	config_lookup_bool(CFG, "psphray.dumpHotspots", &CFG_DUMP_HOTSPOTS);
	config_lookup_bool(CFG, "psphray.traceOnly", &CFG_TRACE_ONLY);
	config_lookup_bool(CFG, "psphray.noPhoton", &CFG_NO_PHOTON);
	config_lookup_bool(CFG, "psphray.HOnly", &CFG_H_ONLY);
	config_lookup_bool(CFG, "psphray.secondaryIonization", &CFG_SECONDARY_IONIZATION);
	config_lookup_bool(CFG, "cosmology.comoving", &CFG_COMOVING);
	config_lookup_bool(CFG, "psphray.disableLTE", &CFG_DISABLE_LTE);
	config_lookup_bool(CFG, "psphray.printRays", &CFG_PRINT_RAYS);

	config_lookup_float(CFG, "psphray.fakeTemperature", &CFG_FAKE_TEMPERATURE);
	config_lookup_float (CFG, "psphray.eos.cloudTemperature", &CFG_EOS_CLOUD_TEMPERATURE);
	config_lookup_float  (CFG, "psphray.recombineThreshold", &CFG_RECOMBINE_THRESHOLD);
	config_lookup_float  (CFG, "psphray.boostSourceFactor", &CFG_BOOST_SOURCE_FACTOR);

	const char * arfilename = NULL;
	const char * xsfilename = NULL;
	const char * eosfilename = NULL;
	config_lookup_string(CFG, "psphray.atomicRates", &arfilename);
	config_lookup_string(CFG, "psphray.xsfilename", &xsfilename);

	ar_init(arfilename);
	xs_init(xsfilename);
	config_lookup_string (CFG, "psphray.eos.filename", &eosfilename);
	if(eosfilename != NULL) {
		eos_init(eosfilename);
		CFG_ENABLE_EOS = CONFIG_TRUE;
	} else {
		CFG_ENABLE_EOS = CONFIG_FALSE;
	}

	lte_init();
	spec_init();
	if(CFG_SECONDARY_IONIZATION)
		secion_init();

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
	if(st != NULL) return st;

	char * pntpath = strdup(path);
	char * dot = rindex(pntpath, '.');
	config_setting_t * pnt = config_root_setting(config);
	char * name = path;
	if(dot != NULL) {
		*dot = 0;
 		name = dot + 1;
		pnt = config_lookup(CFG, pntpath);
		if(pnt == NULL) {
			ERROR("RUNTIME: %s not found", pntpath);
		}
	}
	st = config_setting_add(pnt, name, type);
	free(pntpath);

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
int config_setting_parse_units_elem(config_setting_t * e, int elem, double * rt) {
	config_setting_t * t = config_setting_get_elem(e, elem);
	if(!t) return 0;
	*rt = config_setting_parse_units(t);
	return 1;
}
int config_setting_parse_units_member(config_setting_t * e, char * memb, double * rt) {
	config_setting_t * t = config_setting_get_member(e, memb);
	if(!t) return 0;
	*rt = config_setting_parse_units(t);
	return 1;
}
