
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <messages.h>
#include <array.h>
#include "config.h"
#include "tabfun.h"

int tabfun_setup_col(TabFun * tabfun, char * col, double (*func)(double), double unit) {
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

const double tabfun_get(const TabFun * tabfun, const int id, const double value) {
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

void tabfun_init(TabFun * tabfun, const char * filename) {
	/* the format of the file is stupid fortran. 14 cols per item */
	FILE * fp = fopen(filename, "r");
	if(fp == NULL) {
		ERROR("tabulated function input file %s no access", filename);
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

void tabfun_dump(const TabFun * tabfun, const char * filename) {
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

int tabfun2_ensure_col(TabFun2 * tabfun, char * col, double (*func)(double , double)) {
	int i, j;
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
	*ARRAY_APPEND(tabfun->data, double *) = malloc(sizeof(double) * tabfun->nrows1 * tabfun->nrows2);
	for(i =0; i < tabfun->nrows1; i++) {
		for(j =0; j < tabfun->nrows2; j++) {
		tabfun->data[tabfun->headers_length - 1][i * tabfun->nrows2 + j] = func(tabfun->x[i], tabfun->y[j]);
		}
	}
	return tabfun->headers_length - 1;
}

double tabfun2_get(TabFun2 * tabfun, int id, double x, double y) {
	int xind = (x - tabfun->min1) * tabfun->step1_inv;
	int yind = (y - tabfun->min2) * tabfun->step2_inv;
	int linoffset = -1;
	if(xind < 0) {
		linoffset = 0;
	}
	if(xind >= tabfun->nrows1 - 1) {
		linoffset = (tabfun->nrows1 - 1) * tabfun->nrows2;
	}
	if(linoffset != -1) {
		if(yind < 0) return tabfun->data[id][linoffset + 0];
		if(yind >= tabfun->nrows2 - 1) {
			return tabfun->data[id][linoffset + tabfun->nrows2 - 1];
		}
		double left = tabfun->y[yind];
		double right = tabfun->y[yind+ 1];
		/* the first id is the temprature */
		double leftwt = (right - y);
		double rightwt = (y- left);
		return (leftwt * tabfun->data[id][linoffset + yind] + rightwt * tabfun->data[id][linoffset + yind+1]) * tabfun->step2_inv;
	}
	int lin = -1;
	if(yind < 0) {
		lin = 0;
	}
	if(yind >= tabfun->nrows2 - 1) {
		lin = tabfun->nrows2 - 1;
	}
	if(lin != -1) {
		double left = tabfun->x[xind];
		double right = tabfun->x[xind+ 1];
		/* the first id is the temprature */
		double leftwt = (right - x);
		double rightwt = (x - left);
		return (leftwt * tabfun->data[id][xind * tabfun->nrows2 + lin] + rightwt * tabfun->data[id][(xind + 1) * tabfun->nrows2 + lin]) * tabfun->step1_inv;
	}
	
	/* otherwise safe to do a bilinear */
	double x1 = tabfun->x[xind];
	double x2 = tabfun->x[xind + 1];
	double y1 = tabfun->y[yind];
	double y2 = tabfun->y[yind + 1];
	double fac = tabfun->step1_inv * tabfun->step2_inv;
	double f11 = tabfun->data[id][xind * tabfun->nrows2 + yind];
	double f12 = tabfun->data[id][xind * tabfun->nrows2 + yind + 1];
	double f21 = tabfun->data[id][(xind + 1)* tabfun->nrows2 + yind];
	double f22 = tabfun->data[id][(xind + 1)* tabfun->nrows2 + yind + 1];
	return fac * (f11 * (x2 - x) * (y2 - y) + f21 * (x - x1) * (y2 - y)
			+ f12 * (x2 - x) * (y - y1) + f22 * (x - x1) * (y - y1));
}
