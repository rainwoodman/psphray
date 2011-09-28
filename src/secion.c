#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <messages.h>
#include "config.h"
#include <array.h>

int SECION_PHI_HI;
int SECION_PHI_HEI;
int SECION_EH;

/* Fit taken from 2002 M. Ricotti, N. Y.Gnedin and J.M. Shull ApJ, 575:33-48 */
/* return the normalized values without those prefactors Ej/Ii and such. Aka numbers plotted in figure 6(?) */
/* Also we take the energy in Rydbergs */
/* Notice the thresholds has been reduced by 0.1 to account for discrete erros in the tabulator,
 * making sure at 28 and 11ev the values are not overly artifitially reduced due to the discontinuety */
typedef struct {
	ARRAY_DEFINE_S(data, double *);
	ARRAY_DEFINE_S(headers, char *);
	int nrows1;
	int nrows2;
	int ncols;
	double min1; 
	double min2; 
	double max1; 
	double max2; 
	double step1;
	double step2;
	double step1_inv;
	double step2_inv;
	ARRAY_DEFINE_S(x, double);
	ARRAY_DEFINE_S(y, double);
} TabFun2;

static TabFun2 si = {};

static double a[3][2] = {
	{NAN, 0.2},
	{NAN, 0.2},
	{NAN, 0.4},
};
static double b[3][2] = {
	{0.4092, 0.38},
	{0.4614, 0.38},
	{0.2663, 0.34},
};
static double c[3][2] = {
	{1.7592, 2.0},
	{1.6660, 2.0},
	{1.3163, 2.0},
};
static double C[3][2] = {
	{0.3908, 0.6941},
	{0.0554, 0.0984},
	{1.0, 3.9811},
};

static double secion_get_Phi_factor(int j, double E, double x) {
	if(j >= 2 || j < 0) ERROR("j has to be 0, 1, for HI and HEI respectively");
	double gate = 28 * U_EV / U_RY_ENG;
	int f1 = (E >= (gate-0.1)); /* working out discrete errors in the table */
	double f2 = f1 * pow(gate / E, 0.4);
	
	double y1 = C[j][0] * pow(1. - pow(x, b[j][0]), c[j][0]);
	double y2 = C[j][1] * pow(x, a[j][1]) * pow(1 - pow(x, b[j][1]), c[j][1]);

	return (y1 * f1 - y2 * f2);
}

static double secion_get_Eh_factor(double E, double x) {
	double gate = 11 * U_EV / U_RY_ENG;
	int f1 = (E >= gate - 0.1);
	double f2 = f1 * pow(gate / E, 0.7);
	
	double y1 = C[2][0] * pow(1 - pow(x, b[2][0]), c[2][0]);
	double y2 = C[2][1] * pow(x, a[2][1]) * pow(1 - pow(x, b[2][1]), c[2][1]);

	return (1 - y1 * f1 + y2 * f2);
}

static double Phi_HI(double E, double logx) {
	return secion_get_Phi_factor(0, E, pow(10, logx));
}
static double Phi_HeI(double E, double logx) {
	return secion_get_Phi_factor(1, E, pow(10, logx));
}
static double Eh(double E, double logx) {
	return secion_get_Eh_factor(E, pow(10, logx));
}
static int tabfun2_ensure_col(TabFun2 * tabfun, char * col, double (*func)(double , double));
static double tabfun2_get(TabFun2 * tabfun, int id, double x, double y);

void secion_init() {
	si.nrows1 = 4096;
	si.nrows2 = 1024;
	si.ncols = 3;
	si.min1 = 10 * U_EV / U_RY_ENG;
	si.max1 = 200 * U_EV / U_RY_ENG;
	si.min2 = -8;
	si.max2 = 0;
	si.step1 = (si.max1 - si.min1) / si.nrows1;
	si.step1_inv = 1.0 / si.step1;
	si.step2 = (si.max2 - si.min2) / si.nrows2;
	si.step2_inv = 1.0 / si.step2;
	int i;
	ARRAY_RESIZE(si.x, double, si.nrows1);
	ARRAY_RESIZE(si.y, double, si.nrows2);
	for(i = 0; i < si.nrows1; i++) {
		si.x[i] = si.min1 + i * si.step1;
	}
	for(i = 0; i < si.nrows2; i++) {
		si.y[i] = si.min2 + i * si.step2;
	}
	SECION_PHI_HI = tabfun2_ensure_col(&si, "PhiHI", Phi_HI);
	SECION_PHI_HEI = tabfun2_ensure_col(&si, "PhiHeI", Phi_HeI);
	SECION_EH = tabfun2_ensure_col(&si, "Eh", Eh);
}

double secion_get(int id, double E, double log10x) {
	return tabfun2_get(&si, id, E, log10x);
}

static int tabfun2_ensure_col(TabFun2 * tabfun, char * col, double (*func)(double , double)) {
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
static double tabfun2_get(TabFun2 * tabfun, int id, double x, double y) {
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
