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
	double logTmin; 
	double logTmax; 
	double step;
} AtomicRate;

int AR_LOG_T = -1;
int AR_HI_CI = -1;
int AR_HII_RC_A = -1;
int AR_HII_RCC_A = -1;

AtomicRate ar = {0};


double ar_verner(double Ry) {
	extern float verner_hi_photo_cs_(float *);
	float ry = Ry;
	return verner_hi_photo_cs_(&ry);
}

double ar_get(int id, double logT) {
	if(isnan(logT) || isinf(logT)) 
		ERROR("temperature non-positive: logT = %lg", logT);
	ssize_t index = (logT - ar.logTmin) / ar.step;
	if(index < 0) {
		ERROR("temperature lower than the mininal.(logT= %lg)", logT);
	}
	if(index >= ar.nrows - 1) {
		ERROR("temperature higher than the maximal.(logT= %lg)", logT);
	}
	double left = ar.data[AR_LOG_T][index];
	double right = ar.data[AR_LOG_T][index + 1];
	/* the first col is the temprature */
	double leftwt = (right - logT);
	double rightwt = (logT - left);
	double dlogT = right - left;
	return leftwt * ar.data[id][index] + rightwt * ar.data[id][index+1] / dlogT;
}

static int ar_col_id(char * col) {
	int i;
	for(i = 0; i < ar.ncols; i++) {
		if(!strcasecmp(col, ar.headers[i])) {
			return i;
		}
	}
	ERROR("column %s unknown", col);
}
void ar_init(char * filename) {
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
					&ar.logTmin, &ar.logTmax, &ar.nrows, &ar.step)) {
				ERROR("expecting header logTmin logTmax nrows stepsize, at %s, %d", filename, NR);
			}
			stage++;
		break;
		case 1:
			p = line;
			ar.ncols = len / 14;
			ar.headers = malloc(sizeof(char*) * ar.ncols);
			ar.data = malloc(sizeof(double*) * ar.ncols);
			for(i = 0; i < ar.ncols; i++) {
				p = line + i * 14;
				ar.data[i] = malloc(sizeof(double) * ar.nrows);
				char * q = p + 14 - 1;
				while(*q == ' ' && q >= p) {
					*q = 0;
					q--;
				}
				while(*p == ' ') p++;
				ar.headers[i] = strdup(p);
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
			for(i = 0; i < ar.ncols; i++) {
				p = line + i * 14;
				if(p > line + len) {
					ERROR("file %s line %d is too short", filename, NR);
				}
				p[14 - 1] = 0;
				ar.data[i][irow] = atof(p);
			}
			irow ++;
			if(irow == ar.nrows) stage++;
		break;
		}
		NR++;
		if(stage == 4) break;
	}
	free(line);
	fclose(fp);

	MESSAGE("AR: %d cols, %d rows", ar.ncols, ar.nrows);

	AR_HI_CI = ar_col_id("HIci");
	AR_HII_RC_A = ar_col_id("HIIrcA");
	AR_HII_RCC_A = ar_col_id("HIIrccA");
	AR_LOG_T = 0;
}

void ar_dump(char * filename) {
	FILE * fp = fopen(filename, "w");
	fprintf(fp, "%-.5E %-.5E %-d %-.5E\n", ar.logTmin, ar.logTmax, ar.nrows, ar.step);
	int i;
	for( i = 0; i < ar.ncols; i++) {
		fprintf(fp, " %-14s", ar.headers[i]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	int irow;
	for(irow = 0; irow < ar.nrows; irow++) {
		for( i = 0; i < ar.ncols; i++) {
			fprintf(fp, "%- 14.5E", ar.data[i][irow]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

