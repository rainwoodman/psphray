#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <messages.h>
#include <array.h>

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
int AR_HI_CIC_A = -1;
int AR_HI_CEC_A = -1;
int AR_E_BREMC = -1;
int AR_E_COMPC = -1;
int XS_FREQ = -1;
int XS_HI = -1;

static TabFun ar = {0};
static TabFun xs = {0};

extern float verner_hi_photo_cs_(float *);
static const double tabfun_get(const TabFun * tabfun, const int id, const double value);
static int tabfun_col_id(const TabFun * tabfun, const char * col);
static void tabfun_init(TabFun * tabfun, const char * filename);
static void tabfun_dump(const TabFun * tabfun, const char * filename);
static int tabfun_ensure_col(TabFun * tabfun, char * col, double (*func)(double));

/* analytical forms */
static double verner(double freq) {
	float ryd = freq;
	return verner_hi_photo_cs_(&ryd);
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

	AR_HI_CI = tabfun_ensure_col(&ar, "HIci", NULL);
	AR_HII_RC_A = tabfun_ensure_col(&ar, "HIIrcA", NULL);
	AR_HII_RCC_A = tabfun_ensure_col(&ar, "HIIrccA", NULL);
	AR_HII_RC_B = tabfun_ensure_col(&ar, "HIIrcB", NULL);
	AR_HII_RCC_B = tabfun_ensure_col(&ar, "HIIrccB", NULL);
	AR_HI_CIC_A = tabfun_ensure_col(&ar, "HIcic", NULL);
	AR_HI_CEC_A = tabfun_ensure_col(&ar, "HIcec", NULL);
	AR_E_BREMC = tabfun_ensure_col(&ar, "Brem", bremsstrahlung_cen_1992);
	AR_E_COMPC = tabfun_ensure_col(&ar, "Compton", compton_haiman_1996);

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
		xs.nrows = 1024;
		xs.min = 1;
		xs.max = 16;
		xs.step = (xs.max - xs.min) / (xs.nrows - 1);
		xs.step_inv = 1.0 / xs.step;
		ARRAY_RESIZE(xs.data, double * , xs.ncols);
		ARRAY_RESIZE(xs.headers, char * , xs.ncols);
		xs.headers[0] = "freq";
		xs.data[0] = malloc(sizeof(double)* xs.nrows);
		int i;
		for(i = 0; i < xs.nrows; i++) {
			xs.data[0][i] = xs.min + xs.step * i;
		}
	}

	XS_HI = tabfun_ensure_col(&xs, "HI", verner);
	XS_FREQ = 0;
}

const double xs_get(const int id, const double value) {
	return tabfun_get(&xs, id, value);
}

static int tabfun_ensure_col(TabFun * tabfun, char * col, double (*func)(double)) {
	int i;
	for(i = 0; i < tabfun->headers_length; i++) {
		if(!strcasecmp(col, tabfun->headers[i])) {
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
		tabfun->data[i][j] = func(tabfun->data[0][j]);
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
		index = 0;
//		ERROR("temperature lower than the mininal.(logT= %lg)", logT);
	}
	if(index >= tabfun->nrows - 1) {
		index = tabfun->nrows - 2;
//		ERROR("temperature higher than the maximal.(logT= %lg)", logT);
	}
	float left = tabfun->data[0][index];
	float right = tabfun->data[0][index + 1];
	/* the first col is the temprature */
	float leftwt = (right - value);
	float rightwt = (value- left);
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

