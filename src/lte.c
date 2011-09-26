#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <messages.h>
#include <array.h>
#include <gsl/gsl_randist.h>
#include "config.h"

int LTE_FREQ_HI = 0;
int LTE_FREQ_HEI = 1;
int LTE_FREQ_HEII = 2;

struct {
	ARRAY_DEFINE_S(dist, gsl_ran_discrete_t *);
	size_t Nfreq;
	size_t NlogT;
	double freq_min;
	double freq_max;
	double logT_min;
	double logT_max;
	double logT_step;
	double freq_step;
	ARRAY_DEFINE_S(freq, double);
} lte[3];

static void lte_fill(int id, double freq_th, double freq_max, int xsid);
void lte_init() {
	lte_fill(LTE_FREQ_HI, C_HI_FREQ, 16, XS_HI);
	lte_fill(LTE_FREQ_HEI, C_HEI_FREQ, 16, XS_HEI);
	lte_fill(LTE_FREQ_HEII, C_HEII_FREQ, 16, XS_HEII);
}
double lte_gen_freq(const int id, const double logT) {
	int n = (logT - lte[id].logT_min) / lte[id].logT_step;
	if(n < 0) n = 0;
	if(n >= lte[id].NlogT) n = lte[id].NlogT - 1;
	return lte[id].freq[gsl_ran_discrete(RNG, lte[id].dist[n])];
}

static void lte_fill(int id, double freq_th, double freq_max, int xsid) {
	lte[id].freq_min = freq_th;
	lte[id].freq_max = freq_max;
	lte[id].Nfreq = 1024;
	lte[id].NlogT = 128;
	lte[id].logT_min = 3;
	lte[id].logT_max = 10;
	lte[id].freq_step = (freq_max - freq_th) / lte[id].Nfreq;
	lte[id].logT_step = (lte[id].logT_max - lte[id].logT_min) / lte[id].NlogT;
	ARRAY_RESIZE(lte[id].dist, gsl_ran_discrete_t *, lte[id].NlogT);
	ARRAY_RESIZE(lte[id].freq, double, lte[id].Nfreq);
	int i, j;
	for(j = 0; j < lte[id].Nfreq; j++) {
		double freq = lte[id].freq_min + lte[id].freq_step * j;
		lte[id].freq[j] = freq;
	}
	for(i = 0; i < lte[id].NlogT; i++) {
		double logT = lte[id].logT_min + lte[id].logT_step * i;
		double T = exp(logT);
		double f[lte[id].Nfreq];
		double kT = C_BOLTZMANN * T;
		for(j = 0; j < lte[id].Nfreq; j++) {
			double freq = lte[id].freq[j];
			double sigma = xs_get(xsid, freq);
			f[j] = sigma * freq * freq * freq * exp( - (freq - lte[id].freq_min) * U_RY_ENG / kT);
		}
		lte[id].dist[i] = gsl_ran_discrete_preproc(lte[id].Nfreq, f);
	}
}
