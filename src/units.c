#include <messages.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

double C_OMEGA_L = 0.0;
double C_OMEGA_M = 0.0;
double C_H = 0.0;

double units_simple(char * simple);
double units_init() {
	config_lookup_float(CFG, "cosmology.h", &C_H);
	config_lookup_float(CFG, "cosmology.omegaL", &C_OMEGA_L);
	config_lookup_float(CFG, "cosmology.omegaM", &C_OMEGA_M);
	
	config_lookup_float(CFG, "units.lengthCMh", &U_CM);
	U_CM = C_H / U_CM;

	config_lookup_float(CFG, "units.massGramh", &U_GRAM);
	U_GRAM = C_H / U_GRAM;

	config_lookup_float(CFG, "units.timeSh", &U_SEC);
	U_SEC = C_H / U_SEC;

	U_M = units_simple("m");
	U_KM = units_simple("km");
	U_KPC = units_simple("kpc");
	U_YR = units_simple("yr");
	U_MYR = units_simple("myr");
	U_KG = units_simple("kg");
	U_MSUN = units_simple("msun");
}
double units_simple(char * simple) {
	if(!strcasecmp(simple, "h")) {
		return C_H;
	} else if(!strcasecmp(simple, "cm")) {
		return U_CM;
	} else if(!strcasecmp(simple, "gram")) {
		return U_GRAM;
	} else if(!strcasecmp(simple, "sec")) {
		return U_SEC;
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

double units_format(double value, char * units) {
	double u = units_parse(units);
	return value / u;
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

double z2t(double z) {
	double time;
	double a = 1 / (z + 1);
	const double H0 = 0.1;
	double aeq = pow(C_OMEGA_M / C_OMEGA_L, 1.0/3.0);
	double pre = 2.0 / (3.0 * sqrt(C_OMEGA_L));
	double arg = pow(a / aeq, 3.0 / 2.0)  + sqrt(1.0 + pow(a/ aeq, 3.0));
	return pre * log(arg) / H0;
}
