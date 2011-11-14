#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <messages.h>
#include "config.h"
#include <array.h>
#include "tabfun.h"

int SECION_PHI_HI;
int SECION_PHI_HEI;
int SECION_EH;

/* Fit taken from 2002 M. Ricotti, N. Y.Gnedin and J.M. Shull ApJ, 575:33-48 */
/* return the normalized values without those prefactors Ej/Ii and such. Aka numbers plotted in figure 6(?) */
/* Also we take the energy in Rydbergs */
/* Notice the thresholds has been reduced by 0.1 to account for discrete erros in the tabulator,
 * making sure at 28 and 11ev the values are not overly artifitially reduced due to the discontinuety */
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
	double gate = 28 * U_EV;
	int f1 = (E >= (gate-0.1)); /* working out discrete errors in the table */
	double f2 = f1 * pow(gate / E, 0.4);
	
	double y1 = C[j][0] * pow(1. - pow(x, b[j][0]), c[j][0]);
	double y2 = C[j][1] * pow(x, a[j][1]) * pow(1 - pow(x, b[j][1]), c[j][1]);

	return (y1 * f1 - y2 * f2);
}

static double secion_get_Eh_factor(double E, double x) {
	double gate = 11 * U_EV;
	int f1 = (E >= gate - 0.1);
	double f2 = f1 * pow(gate / E, 0.7);
	
	double y1 = C[2][0] * pow(1 - pow(x, b[2][0]), c[2][0]);
	double y2 = C[2][1] * pow(x, a[2][1]) * pow(1 - pow(x, b[2][1]), c[2][1]);

	return (1. - y1 * f1 + y2 * f2);
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

void secion_init() {
	si.nrows1 = 4096;
	si.nrows2 = 1024;
	si.ncols = 3;
	si.min1 = 10 * U_EV;
	si.max1 = 200 * U_EV;
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
	double rt = tabfun2_get(&si, id, E, log10x);
	return rt;
}

