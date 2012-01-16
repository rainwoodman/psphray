#include <messages.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"
double U_CM = 0.0;
double U_GRAM = 0.0;
double U_SEC = 0.0;
double U_HZ = 0.0;
double U_M = 0.0;
double U_KM = 0.0;
double U_KPC = 0.0;
double U_YR = 0.0;
double U_MYR = 0.0;
double U_KG = 0.0;
double U_MSUN = 0.0;
double U_MPROTON = 0.0;
double U_KELVIN = 0.0;

double U_J = 0.0;
double U_ERG = 0.0;
double U_EV = 0.0;
double U_KEV = 0.0;
double U_RY_ENG = 0.0;

double C_BOLTZMANN = 0.0;
double C_PLANCK = 0.0;
double C_SPEED_LIGHT = 0.0;

double C_H_PER_MASS = 0.0;
double C_HE_PER_MASS = 0.0;

double C_HI_ENERGY = 0.0;
double C_HEI_ENERGY = 0.0;
double C_HEII_ENERGY = 0.0;
double C_SOLAR_LUM = 0.0;

double C_GRAVITY = 0.0;
double C_HUBBLE = 0.0;
double U_MPROTON_OVER_C_BOLTZMANN = 0.0;
double C_BOLTZMANN_OVER_U_MPROTON = 0.0;

	static struct {
		const char * name;
		double * value;
	} u_list[] = {
		{"h", &C_H},
		{"hp", &C_PLANCK},
		{"kb", &C_BOLTZMANN},
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
		{"ev", &U_EV},
		{"kev", &U_KEV},
		{NULL, NULL},
	};
double units_simple(char * simple);
double units_init() {

	U_CM = C_H / C_1_CMH;

	U_GRAM = C_H / C_1_GRAMH;

	U_SEC = C_H / C_1_SECH;

	U_HZ = 1.0 / U_SEC;
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
	C_PLANCK = 6.626068e-34 * U_M * U_M * U_KG / U_SEC;
	C_GRAVITY = 6.674e-11 * U_M * U_M * U_M / U_KG / U_SEC / U_SEC;
	C_HUBBLE = C_H * 100 * 1e3 * U_M / U_SEC / (1000 * U_KPC);
	U_J = U_KG * U_M * U_M / (U_SEC * U_SEC);
	U_EV = 1.602176487e-19 * U_J;
	U_KEV = U_EV * 1000;
	U_ERG = 1e-7 * U_J;
	U_RY_ENG = 13.60569253 * U_EV;
	C_H_PER_MASS = C_HMF / U_MPROTON;
	C_HE_PER_MASS = C_HEMF / (4 * U_MPROTON);
	C_SPEED_LIGHT = 3e8 * U_M / U_SEC;

	C_HI_ENERGY = 13.60569253 * U_EV;
	C_HEI_ENERGY = 24.587 * U_EV;
	C_HEII_ENERGY = 54.416 * U_EV;

	C_SOLAR_LUM = 3.899e33 * U_ERG / U_SEC;
	U_MPROTON_OVER_C_BOLTZMANN = U_MPROTON / C_BOLTZMANN;
	C_BOLTZMANN_OVER_U_MPROTON = 1.0 / (U_MPROTON / C_BOLTZMANN);
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
double t2z(double t) {
	double a;
	const double H0 = 0.1;
	double aeq = pow(C_OMEGA_M / C_OMEGA_L, 1.0/3.0);
	double pre = 2.0 / (3.0 * sqrt(C_OMEGA_L));
	a = pow(sinh(H0 * t / pre), 2.0/3.0) * aeq;
	return 1.0 / a - 1.0;
}

