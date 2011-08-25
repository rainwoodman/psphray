#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <messages.h>

typedef struct {
	double ** data;
	char ** headers;
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

int XS_FREQ = -1;
int XS_HI = -1;

static TabFun ar = {0};
static TabFun xs = {0};

extern float verner_hi_photo_cs_(float *);
static const double tabfun_get(const TabFun * tabfun, const int id, const double value);
static int tabfun_col_id(const TabFun * tabfun, const char * col);
static void tabfun_init(TabFun * tabfun, const char * filename);
static void tabfun_dump(const TabFun * tabfun, const char * filename);

void ar_init(const char * filename) {
	tabfun_init(&ar, filename);
	MESSAGE("AR: %d cols, %d rows", ar.ncols, ar.nrows);

	AR_HI_CI = tabfun_col_id(&ar, "HIci");
	AR_HII_RC_A = tabfun_col_id(&ar, "HIIrcA");
	AR_HII_RCC_A = tabfun_col_id(&ar, "HIIrccA");
	AR_HII_RC_B = tabfun_col_id(&ar, "HIIrcB");
	AR_HII_RCC_B = tabfun_col_id(&ar, "HIIrccB");
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
		xs.ncols = 2;
		xs.nrows = 1024;
		xs.min = 1;
		xs.max = 16;
		xs.step = (xs.max - xs.min) / (xs.nrows - 1);
		xs.step_inv = 1.0 / xs.step;
		xs.data = malloc(sizeof(double*)* xs.ncols);
		xs.headers = malloc(sizeof(char*)* xs.ncols);
		xs.headers[0] = "freq";
		xs.headers[1] = "HI";
		xs.data[0] = malloc(sizeof(double)* xs.nrows);
		xs.data[1] = malloc(sizeof(double)* xs.nrows);
		int i;
		for(i = 0; i < xs.nrows; i++) {
			xs.data[0][i] = xs.min + xs.step * i;
			float ry = xs.data[0][i];
			xs.data[1][i] = verner_hi_photo_cs_(&ry);
		}
	}

	XS_HI = tabfun_col_id(&xs, "HI");
	XS_FREQ = 0;
}

const double xs_get(const int id, const double value) {
	return tabfun_get(&xs, id, value);
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

static int tabfun_col_id(const TabFun * tabfun, const char * col) {
	int i;
	for(i = 0; i < tabfun->ncols; i++) {
		if(!strcasecmp(col, tabfun->headers[i])) {
			return i;
		}
	}
	ERROR("column %s unknown", col);
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
			tabfun->headers = malloc(sizeof(char*) * tabfun->ncols);
			tabfun->data = malloc(sizeof(double*) * tabfun->ncols);
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
	for( i = 0; i < tabfun->ncols; i++) {
		fprintf(fp, " %-14s", tabfun->headers[i]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	int irow;
	for(irow = 0; irow < tabfun->nrows; irow++) {
		for( i = 0; i < tabfun->ncols; i++) {
			fprintf(fp, "%- 14.5E", tabfun->data[i][irow]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

