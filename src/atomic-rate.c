#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <messages.h>
#include <array.h>
#include "config.h"

#include "tabfun.h"

int AR_LOG_T = -1;
int AR_HI_CI = -1;
int AR_HII_RC_A = -1;
int AR_HII_RCC_A = -1;
int AR_HII_RC_B = -1;
int AR_HII_RCC_B = -1;
int AR_HI_CIC = -1;
int AR_HI_CEC = -1;
int AR_E_BREMC = -1;
int AR_E_COMPC = -1;
int AR_HEI_CI = -1;
int AR_HEII_CI = -1;
int AR_HEII_RC_A = -1;
int AR_HEII_RC_B = -1;
int AR_HEIII_RC_A = -1;
int AR_HEIII_RC_B = -1;
int XS_ENG = -1;
int XS_HI = -1;
int XS_HEI = -1;
int XS_HEII = -1;

static TabFun ar = {0};
static TabFun xs = {0};
/* holding the low cuts of the xsections. the analytical form doesn't have this cut, because
 * tabulated function can't do the sharp cut well. 100 here is just a number greater than
 * the max possible number of crosssctions (which is 3 for now)*/
static double xs_cut[100];


extern float verner_hi_photo_cs_(float *);
extern float osterbrok_hei_photo_cs_(float *);
extern float osterbrok_heii_photo_cs_(float *);

/* analytical forms */
static double verner(double eng) {
	float ryd = eng / C_HI_ENERGY;
	return verner_hi_photo_cs_(&ryd);
}
static double osterbrok_HeI(double eng) {
	float ryd = eng / C_HI_ENERGY;
	return osterbrok_hei_photo_cs_(&ryd);
}
static double osterbrok_HeII(double eng) {
	float ryd = eng / C_HI_ENERGY;
	return osterbrok_heii_photo_cs_(&ryd);
}
static double bremsstrahlung_cen_1992(double T) {
	return 1.42e-27 * 1.5 * sqrt(T);
}
/* not there yet, need to find a way to deal with the backgroun temperature */
static double compton_haiman_1996(double T) {
	return 0.0;
}

void ar_init(const char * filename) {
	tabfun_init(&ar, filename);
	MESSAGE("AR: %d cols, %d rows", ar.ncols, ar.nrows);

	AR_HI_CI = tabfun_setup_col(&ar, "HIci", NULL, U_CM * U_CM * U_CM / U_SEC);
	AR_HII_RC_A = tabfun_setup_col(&ar, "HIIrcA", NULL, U_CM * U_CM * U_CM / U_SEC);
	AR_HII_RCC_A = tabfun_setup_col(&ar, "HIIrccA", NULL, U_ERG * U_CM * U_CM * U_CM / U_SEC);
	AR_HII_RC_B = tabfun_setup_col(&ar, "HIIrcB", NULL, U_CM * U_CM * U_CM / U_SEC);
	AR_HII_RCC_B = tabfun_setup_col(&ar, "HIIrccB", NULL, U_ERG * U_CM * U_CM * U_CM / U_SEC);
	AR_HEI_CI = tabfun_setup_col(&ar, "HeIci", NULL, U_CM * U_CM * U_CM / U_SEC);
	AR_HEII_CI = tabfun_setup_col(&ar, "HeIIci", NULL, U_CM * U_CM * U_CM / U_SEC);
	AR_HEII_RC_A = tabfun_setup_col(&ar, "HeIIrcA", NULL, U_CM * U_CM * U_CM / U_SEC);
	AR_HEII_RC_B = tabfun_setup_col(&ar, "HeIIrcB", NULL, U_CM * U_CM * U_CM / U_SEC);
	AR_HEIII_RC_A = tabfun_setup_col(&ar, "HeIIIrcA", NULL, U_CM * U_CM * U_CM / U_SEC);
	AR_HEIII_RC_B = tabfun_setup_col(&ar, "HeIIIrcB", NULL, U_CM * U_CM * U_CM / U_SEC);
	AR_HI_CIC = tabfun_setup_col(&ar, "HIcic", NULL, U_ERG * (U_CM * U_CM * U_CM) / U_SEC);
	AR_HI_CEC = tabfun_setup_col(&ar, "HIcec", NULL, U_ERG * (U_CM * U_CM * U_CM) / U_SEC);
	AR_E_BREMC = tabfun_setup_col(&ar, "Brem", bremsstrahlung_cen_1992, U_ERG * (U_CM * U_CM * U_CM) / U_SEC);
	AR_E_COMPC = tabfun_setup_col(&ar, "Compton", compton_haiman_1996, U_ERG / U_SEC);

	AR_LOG_T = 0;
}

const double ar_get(const int id, const double value) {
	return tabfun_get(&ar, id, value);
}

void xs_init(const char * filename) {
	if(filename != NULL) {
		tabfun_init(&xs, filename);
		MESSAGE("XS: %d cols, %d rows", xs.ncols, xs.nrows);
	} else {
		xs.ncols = 1;
		xs.nrows = 16384;
		xs.min = 1 * C_HI_ENERGY;
		xs.max = 100 * C_HI_ENERGY;
		xs.step = (xs.max - xs.min) / (xs.nrows - 1);
		xs.step_inv = 1.0 / xs.step;
		ARRAY_RESIZE(xs.data, double * , xs.ncols);
		ARRAY_RESIZE(xs.headers, char * , xs.ncols);
		xs.headers[0] = "eng";
		xs.data[0] = malloc(sizeof(double)* xs.nrows);
		int i;
		for(i = 0; i < xs.nrows; i++) {
			xs.data[0][i] = xs.min + xs.step * i;
		}
	}

	XS_HI = tabfun_setup_col(&xs, "HI", verner, U_CM * U_CM);
	XS_HEI = tabfun_setup_col(&xs, "HeI", osterbrok_HeI, U_CM * U_CM);
	XS_HEII = tabfun_setup_col(&xs, "HeII", osterbrok_HeII, U_CM * U_CM);
	/* deliberately not use the exact value because the input is
     * less exact sometimes is ugly. aka we don't want to ionizize nothing
     * when the spectra is a 13.6 ev mono. */
	xs_cut[XS_HI] = 13.6 * U_EV; //C_HI_ENERGY;  
	xs_cut[XS_HEI] = 24.5 * U_EV; //C_HEI_ENERGY;
	xs_cut[XS_HEII] = 54.4 * U_EV; //C_HEII_ENERGY;
	XS_ENG = 0;
}

const double xs_get(const int id, const double value) {
	if(value < xs_cut[id]) {
		return 0.0;
	}
	return tabfun_get(&xs, id, value);
}

