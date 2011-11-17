#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <messages.h>
#include <array.h>
#include "config.h"

#include "tabfun.h"

static TabFun eos = {0};

static int EOS_X;
static int EOS_EGYHOT;
void eos_init(const char * filename) {
	tabfun_init(&eos, filename);
	MESSAGE("AR: %d cols, %d rows", eos.ncols, eos.nrows);

	EOS_X = tabfun_setup_col(&eos, "x", NULL, 1);
	EOS_EGYHOT = tabfun_setup_col(&eos, "egyhot", NULL, 1);
}

double eos_get_cloud_fraction(const double dens_phys) {
	double ll = log10(dens_phys);
	if(ll <= eos.data[0][0]) return 0.0;
	return tabfun_get(&eos, EOS_X, ll);
}
double eos_get_egyhot(const double dens_phys) {
	double ll = log10(dens_phys);
	if(ll <= eos.data[0][0]) return 0.0;
	return tabfun_get(&eos, EOS_EGYHOT, ll);
}
