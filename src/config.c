#include <stdio.h>
#include <string.h>
#include <messages.h>
#include <libconfig/libconfig.h>

extern void units_init();

config_t CFG[] = {{0}};

config_setting_t * config_ensure(config_t * config, char * path, int type);
#define config_ensure_int(c, p, v)  if(!config_lookup(c, p)) config_setting_set_int(config_ensure(c, p, CONFIG_TYPE_INT), v)
#define config_ensure_int64(c, p, v)  if(!config_lookup(c,p)) config_setting_set_int64(config_ensure(c, p, CONFIG_TYPE_INT64), v)
#define config_ensure_float(c, p, v)  if(!config_lookup(c,p)) config_setting_set_float(config_ensure(c, p, CONFIG_TYPE_FLOAT), v)
#define config_ensure_bool(c, p, v)  if(!config_lookup(c,p)) config_setting_set_bool(config_ensure(c, p, CONFIG_TYPE_BOOL), v)
#define config_ensure_string(c, p, v)  if(!config_lookup(c,p)) config_setting_set_string(config_ensure(c, p, CONFIG_TYPE_STRING), v)

void cfg_init(char * filename) {
	config_init(CFG);
	if(!config_read_file(CFG, filename)) {
		ERROR("%s: %d: %s", config_error_file(CFG), config_error_line(CFG), config_error_text(CFG));
	}
	config_ensure        (CFG, "psphray", CONFIG_TYPE_GROUP);
	config_ensure_string (CFG, "psphray.atomicRates", "<filename>");
	config_ensure        (CFG, "cosmology", CONFIG_TYPE_GROUP);
	config_ensure_bool   (CFG, "cosmology.comoving", CONFIG_TRUE);
	config_ensure_float  (CFG, "cosmology.omegaL", 0.74);
	config_ensure_float  (CFG, "cosmology.omegaM", 0.26);
	config_ensure_float  (CFG, "cosmology.omegaB", 0.044);
	config_ensure_float  (CFG, "cosmology.h", 0.72);
	config_ensure        (CFG, "box", CONFIG_TYPE_GROUP);
	config_ensure_float  (CFG, "box.boxsize", 0.0);
	config_ensure_string (CFG, "box.boundary", "vaccum");
	config_ensure        (CFG, "units", CONFIG_TYPE_GROUP);
	config_ensure_float  (CFG, "units.massGramh", 1.9847005219450705e+43);
	config_ensure_float  (CFG, "units.lengthCMh", 3.0835455558318480e+21);
	config_ensure_float  (CFG, "units.timeSh", 3.08568025e+16);
	config_ensure        (CFG, "outputs", CONFIG_TYPE_GROUP);
	config_ensure_string (CFG, "outputs.mode", "fixed");
	config_ensure_string (CFG, "outputs.prefix", "../test/");
	config_ensure        (CFG, "outputs.fixed", CONFIG_TYPE_GROUP);
	config_ensure_int    (CFG, "outputs.fixed.nsteps", 5);
	config_ensure        (CFG, "epochs", CONFIG_TYPE_LIST);
	units_init();
}
void cfg_dump(char * filename) {
	config_write_file(CFG, filename);
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
