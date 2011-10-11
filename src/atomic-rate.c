#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <messages.h>
#include <array.h>
#include "config.h"

typedef struct {
	ARRAY_DEFINE_S(data, double *);
	ARRAY_DEFINE_S(headers, char *);
	int nrows;
	int ncols;
	double min; 
	double max; 
	double step;
	double step_inv;
} TabFun;

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

static const double tabfun_get(const TabFun * tabfun, const int id, const double value);
static int tabfun_col_id(const TabFun * tabfun, const char * col);
static void tabfun_init(TabFun * tabfun, const char * filename);
static void tabfun_dump(const TabFun * tabfun, const char * filename);
static int tabfun_setup_col(TabFun * tabfun, char * col, double (*func)(double), double unit);

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

static int tabfun_setup_col(TabFun * tabfun, char * col, double (*func)(double), double unit) {
/* if the col exists, scale them by the unit. if not, calculate from the given function func and scale by unit */
	int i;
	for(i = 0; i < tabfun->headers_length; i++) {
		if(!strcasecmp(col, tabfun->headers[i])) {
			int j;
			for(j =0; j < tabfun->nrows; j++) {
				tabfun->data[i][j] *= unit;
			}
			return i;
		}
	}
	if(func == NULL) {
		ERROR("column %s unknown and no hard coded analytical form", col);
	} else {
		MESSAGE("using analytical form of column %s", col);
	}
	*ARRAY_APPEND(tabfun->headers, char*) = strdup(col);
	*ARRAY_APPEND(tabfun->data, double *) = malloc(sizeof(double) * tabfun->nrows);
	int j;
	for(j =0; j < tabfun->nrows; j++) {
		tabfun->data[i][j] = func(tabfun->data[0][j]) * unit;
	}
	return i;
}

static const double tabfun_get(const TabFun * tabfun, const int id, const double value) {
/*
	if(isnan(logT) || isinf(logT)) 
		ERROR("temperature non-positive: logT = %lg", logT);
*/
	ssize_t index = (value - tabfun->min) * tabfun->step_inv;
	if(index < 0) {
//		ERROR("temperature lower than the mininal.(logT= %lg)", logT);
		return tabfun->data[id][0];
	}
	if(index >= tabfun->nrows - 1) {
		return tabfun->data[id][tabfun->nrows - 1];
//		ERROR("temperature higher than the maximal.(logT= %lg)", logT);
	}
	double left = tabfun->data[0][index];
	double right = tabfun->data[0][index + 1];
	/* the first col is the temprature */
	double leftwt = (right - value);
	double rightwt = (value- left);
	return (leftwt * tabfun->data[id][index] + rightwt * tabfun->data[id][index+1]) * tabfun->step_inv;
}

static void tabfun_init(TabFun * tabfun, const char * filename) {
	FILE * fp = fopen(filename, "r");
	if(fp == NULL) {
		ERROR("atomic rates %s no access", filename);
	}
	int NR = 0;
	char * line = NULL;
	size_t size = 0;
	size_t len;
	int stage = 0;
	int irow = 0;
	char * p;
	int i;
	while(0 <= (len = getline(&line, &size, fp))) {
		if(line[0] == '#') {
			NR++;
			continue;
		}
		switch(stage) {
		case 0:
			if(4 != sscanf(line, "%lf %lf %d %lf", 
					&tabfun->min, &tabfun->max, &tabfun->nrows, &tabfun->step)) {
				ERROR("expecting header logTmin logTmax nrows stepsize, at %s, %d", filename, NR);
			}
			tabfun->step_inv = 1.0 / tabfun->step;
			stage++;
		break;
		case 1:
			p = line;
			tabfun->ncols = len / 14;
			ARRAY_RESIZE(tabfun->data, double * , tabfun->ncols);
			ARRAY_RESIZE(tabfun->headers, char * , tabfun->ncols);
			for(i = 0; i < tabfun->ncols; i++) {
				p = line + i * 14;
				tabfun->data[i] = malloc(sizeof(double) * tabfun->nrows);
				char * q = p + 14 - 1;
				while(*q == ' ' && q >= p) {
					*q = 0;
					q--;
				}
				while(*p == ' ') p++;
				tabfun->headers[i] = strdup(p);
			}
			stage++;
		break;
		case 2:
			/* skip the next line citing the source of the data */
			stage++;
			irow = 0;
		break;
		case 3:
			p = line;
			for(i = 0; i < tabfun->ncols; i++) {
				p = line + i * 14;
				if(p > line + len) {
					ERROR("file %s line %d is too short", filename, NR);
				}
				p[14 - 1] = 0;
				tabfun->data[i][irow] = atof(p);
			}
			irow ++;
			if(irow == tabfun->nrows) stage++;
		break;
		}
		NR++;
		if(stage == 4) break;
	}
	free(line);
	fclose(fp);
}

static void tabfun_dump(const TabFun * tabfun, const char * filename) {
	FILE * fp = fopen(filename, "w");
	fprintf(fp, "%-.5E %-.5E %-d %-.5E\n", tabfun->min, tabfun->max, tabfun->nrows, tabfun->step);
	int i;
	for( i = 0; i < tabfun->headers_length; i++) {
		fprintf(fp, " %-14s", tabfun->headers[i]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	int irow;
	for(irow = 0; irow < tabfun->nrows; irow++) {
		for( i = 0; i < tabfun->headers_length; i++) {
			fprintf(fp, "%- 14.5E", tabfun->data[i][irow]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

