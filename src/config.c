#include <stdio.h>
#include <string.h>
#include <messages.h>
#include <libconfig/libconfig.h>

static config_t cfg[] = {{0}};

config_setting_t * config_ensure(config_t * config, char * path, int type);
#define config_ensure_int(c, p, v)  if(!config_lookup(c, p)) config_setting_set_int(config_ensure(c, p, CONFIG_TYPE_INT), v)
#define config_ensure_int64(c, p, v)  if(!config_lookup(c,p)) config_setting_set_int64(config_ensure(c, p, CONFIG_TYPE_INT64), v)
#define config_ensure_float(c, p, v)  if(!config_lookup(c,p)) config_setting_set_float(config_ensure(c, p, CONFIG_TYPE_FLOAT), v)
#define config_ensure_bool(c, p, v)  if(!config_lookup(c,p)) config_setting_set_bool(config_ensure(c, p, CONFIG_TYPE_BOOL), v)
#define config_ensure_string(c, p, v)  if(!config_lookup(c,p)) config_setting_set_string(config_ensure(c, p, CONFIG_TYPE_STRING), v)

void cfg_init(char * filename) {
	config_init(cfg);
	if(!config_read_file(cfg, filename)) {
		ERROR("%s: %d: %s", config_error_file(cfg), config_error_line(cfg), config_error_text(cfg));
	}
	config_ensure        (cfg, "psphray", CONFIG_TYPE_GROUP);
	config_ensure_string (cfg, "psphray.atomicRates", "<filename>");
	config_ensure        (cfg, "cosmology", CONFIG_TYPE_GROUP);
	config_ensure_bool   (cfg, "cosmology.comoving", CONFIG_TRUE);
	config_ensure_float  (cfg, "cosmology.omegaL", 0.74);
	config_ensure_float  (cfg, "cosmology.omegaM", 0.26);
	config_ensure_float  (cfg, "cosmology.omegaB", 0.044);
	config_ensure_float  (cfg, "cosmology.h", 0.72);
	config_ensure        (cfg, "box", CONFIG_TYPE_GROUP);
	config_ensure_float  (cfg, "box.boxsize", 0.0);
	config_ensure_string (cfg, "box.boundary", "vaccum");
	config_ensure        (cfg, "units", CONFIG_TYPE_GROUP);
	config_ensure_float  (cfg, "units.massGramh", 1.9847005219450705e+43);
	config_ensure_float  (cfg, "units.lengthCMh", 3.0835455558318480e+21);
	config_ensure_float  (cfg, "units.timeSh", 3.08568025e+16);
	config_ensure        (cfg, "outputs", CONFIG_TYPE_GROUP);
	config_ensure_string (cfg, "outputs.mode", "fixed");
	config_ensure_string (cfg, "outputs.prefix", "../test/");
	config_ensure        (cfg, "outputs.fixed", CONFIG_TYPE_GROUP);
	config_ensure_int    (cfg, "outputs.fixed.nsteps", 5);
	config_ensure        (cfg, "epochs", CONFIG_TYPE_LIST);
}
void cfg_dump(char * filename) {
	config_write_file(cfg, filename);
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
			pnt = config_lookup(cfg, pntpath);
			if(pnt == NULL) {
				ERROR("RUNTIME: %s not found", pntpath);
			}
			st = config_setting_add(pnt, name, type);
		}
		free(pntpath);
	}
	return st;
}
