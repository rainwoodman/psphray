#include <uthash/uthash.h>
#include <messages.h>
#include "config.h"

double U_CM = 0.0;
double U_GRAM = 0.0;
double U_SEC = 0.0;
double U_M = 0.0;
double U_KM = 0.0;
double U_KPC = 0.0;
double U_YR = 0.0;
double U_MYR = 0.0;
double U_KG = 0.0;
double U_MSUN = 0.0;

double units_simple(char * simple);
double units_init() {
	U_CM = units_simple("cm");
	U_GRAM = units_simple("gram");
	U_SEC = units_simple("sec");
	U_M = units_simple("m");
	U_KM = units_simple("km");
	U_KPC = units_simple("kpc");
	U_YR = units_simple("yr");
	U_MYR = units_simple("myr");
	U_KG = units_simple("kg");
	U_MSUN = units_simple("msun");

}
double units_simple(char * simple) {
	double h;
	config_lookup_float(CFG, "cosmology.h", &h);
	if(!strcasecmp(simple, "h")) {
		return h;
	} else if(!strcasecmp(simple, "cm")) {
		double t;
		config_lookup_float(CFG, "units.lengthCMh", &t);
		return h / t;
	} else if(!strcasecmp(simple, "gram")) {
		double t;
		config_lookup_float(CFG, "units.massGramh", &t);
		return h / t;
	} else if(!strcasecmp(simple, "sec")) {
		double t;
		config_lookup_float(CFG, "units.timeSh", &t);
		return h / t;
	} else if(!strcasecmp(simple, "m")) {
		return units_simple("cm") * 100 ;
	} else if(!strcasecmp(simple, "km")) {
		return units_simple("m") * 1000 ;
	} else if(!strcasecmp(simple, "kpc")) {
		return units_simple("m") * 3.08568025e19;
	} else if(!strcasecmp(simple, "yr")) {
		return units_simple("sec") * 31556926.0;
	} else if(!strcasecmp(simple, "kg")) {
		return units_simple("gram") * 1000 ;
	} else if(!strcasecmp(simple, "myr")) {
		return units_simple("yr") * 1e6 ;
	} else if(!strcasecmp(simple, "msun")) {
		return units_simple("kg") * 1.98892e30;
	}
	ERROR("unit %s unknown", simple);
}

double units_parse(char * units) {
	units = strdup(units);
	char * p, *token;
	char op = '*';
	token = units;
	p = units;
	double unit = 1.0;
	double section = 1.0;
	int skip_blank = 1;
	while(*p) {
		switch(*p) {
		case ' ':
			if(skip_blank) {
				p++;
				continue;
			}
			char save = *p;
			*p = 0;
			section *= units_simple(token);
			*p = save;
			skip_blank = 1;
			p++;
		break;
		case '/':
		case '*':
			switch(op) {
				case '/': unit /= section; break;
				case '*': unit *= section; break;
			}
			section = 1.0;
			op = *p;
			p++;
			skip_blank = 1;
		break;
		default:
			if(skip_blank) {
				token = p;
				skip_blank = 0;
			}
			p++;
		}
	}
	section *= units_simple(token);
	switch(op) {
		case '/': unit /= section; break;
		case '*': unit *= section; break;
	}
	free(units);
	return unit;
}

