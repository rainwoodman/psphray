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
double U_MPROTON = 0.0;
double U_KELVIN = 0.0;

double C_OMEGA_L = 0.0;
double C_OMEGA_M = 0.0;
double C_H = 0.0;
double C_HMF = 0.0;
double C_BOLTZMANN = 0.0;

	static struct {
		const char * name;
		double * value;
	} u_list[] = {
		{"h", &C_H},
		{"cm", &U_CM},
		{"m", &U_M},
		{"km", &U_KM},
		{"kpc", &U_KPC},
		{"g", &U_GRAM},
		{"kg", &U_KG},
		{"yr", &U_YR},
		{"myr", &U_MYR},
		{"msun", &U_MSUN},
		{"mproton", &U_MPROTON},
		{"kelvin", &U_KELVIN},
		{NULL, NULL},
	};
double units_simple(char * simple);
double units_init() {
	config_lookup_float(CFG, "cosmology.h", &C_H);
	config_lookup_float(CFG, "cosmology.hmf", &C_HMF);
	config_lookup_float(CFG, "cosmology.omegaL", &C_OMEGA_L);
	config_lookup_float(CFG, "cosmology.omegaM", &C_OMEGA_M);
	
	config_lookup_float(CFG, "units.lengthCMh", &U_CM);
	U_CM = C_H / U_CM;

	config_lookup_float(CFG, "units.massGramh", &U_GRAM);
	U_GRAM = C_H / U_GRAM;

	config_lookup_float(CFG, "units.timeSh", &U_SEC);
	U_SEC = C_H / U_SEC;

	U_M = U_CM * 100;
	U_KM = U_M * 1000;
	U_KPC = U_M * 3.08568025e19;
	U_YR = U_SEC * 31556926.0;
	U_KG = U_GRAM * 1000;
	U_MYR = U_YR * 1e6;
	U_MSUN = U_KG * 1.98892e30;
	U_MPROTON = U_KG * 1.67262158e-27;
	U_KELVIN = 1.0;
	C_BOLTZMANN = 1.3806503e-23 * U_M * U_M * U_KG / U_SEC / U_SEC;
}

double units_simple(char * simple) {
	int i;
	for (i = 0; i < sizeof(u_list) / sizeof(u_list[0]); i++) {
		if(!strcasecmp(simple, u_list[i].name)) return * u_list[i].value;
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

float ieye2T(const float ie, const float ye) {
	return U_MPROTON / C_BOLTZMANN * ie / ( (1.0 - C_HMF) * 0.25 + 
			C_HMF + C_HMF * ye) * 2.0 / 3.0 ;
}
