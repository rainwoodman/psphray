#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <bitmask.h>
#include <messages.h>

#include "config.h"
#include "psystem.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

struct x_t {
	intptr_t ipar;
	float d;
	float b;
};

struct r_t {
	float s[3];
	float dir[3];
	double Nph;
	double freq;
	intptr_t * ipars;
	struct x_t * x;
	double length;
	size_t ipars_size;
	size_t ipars_length;
	intptr_t isrc;
};

extern PSystem psys;
typedef struct _Solver Solver;
extern Solver * solver_new();
extern int solver_evolve(Solver * s, intptr_t ipar);
extern void solver_delete(Solver * s);
extern float sph_depth(const float r_h);

static int x_t_compare(const void * p1, const void * p2);
static float dist(const float p1[3], const float p2[3]);

static void evolve_srcs(struct r_t * r, size_t Nr);
static void trace(struct r_t * r, size_t Nr);
static size_t deposit(struct r_t * r, size_t Nr, intptr_t **ipars, size_t * ipars_size);
static void solve(intptr_t * ipars, size_t ipars_length);
static size_t recombine_rays(intptr_t * ipars, size_t ipars_length, struct r_t ** r, size_t * r_size, size_t extra);


gsl_rng * rng = NULL;;
double * reservoir = NULL;
char * active = NULL;
void init() {
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, 123456);
}

void run() {
	size_t Nrecrays = 0;
	size_t Nsrcrays = 1;
	struct r_t * r = NULL;
	size_t r_size = 0;
	intptr_t * ipars = NULL;
	size_t ipars_size = 0;	
	size_t ipars_length = 0;
	/* */

	active = bitmask_alloc(psys.npar);
	bitmask_clear_all(active);
	reservoir = calloc(sizeof(double), psys.nsrcs);

	while(psys.tick < psys.nticks) {
		Nrecrays = recombine_rays(ipars, ipars_length, &r, &r_size, Nsrcrays);
		evolve_srcs(r, Nsrcrays);
		MESSAGE("rays: %lu(rec) %lu(src)", Nrecrays, Nsrcrays);

		/* trace rays */
		trace(r, Nrecrays + Nsrcrays);
		/* deposit photons */
		ipars_length = deposit(r, Nrecrays + Nsrcrays, &ipars, &ipars_length);
		MESSAGE("pars: %lu", ipars_length);

		solve(ipars, ipars_length);

		psys.tick++;
	}

	free(ipars);
	free(r);
	free(reservoir);
	free(active);
}

static void evolve_srcs(struct r_t * r, size_t Nr) {
	intptr_t i;
	
	gsl_ran_discrete_t * randist = gsl_ran_discrete_preproc(psys.nsrcs, reservoir);

	size_t raycount[psys.nsrcs];

	for(i = 0; i < psys.nsrcs; i++) {
		raycount[i] = 0;
		reservoir[i] += psys.srcs[i].Ngamma_sec * psys.tick_time / U_SEC;
	}
	
	for(i = 0; i < Nr; i++) {
		intptr_t isrc = gsl_ran_discrete(rng, randist);
		r[i].isrc = isrc;
		raycount[isrc]++;
	}
	for(i = 0; i < Nr; i++) {
		intptr_t isrc = r[i].isrc;
		int d;
		for(d = 0; d < 3; d++) {
			r[i].s[d] = psys.srcs[isrc].pos[d];
		}
		double dx, dy, dz;
		gsl_ran_dir_3d(rng, &dx, &dy, &dz);
		r[i].dir[0] = dx;
		r[i].dir[1] = dy;
		r[i].dir[2] = dz;
		r[i].Nph = reservoir[isrc] / raycount[isrc];
		r[i].freq = 1.0;
		r[i].length = psys.boxsize;
	}
	for(i = 0; i < Nr; i++) {
		intptr_t isrc = r[i].isrc;
		reservoir[isrc] = 0.0;
	}
	gsl_ran_discrete_free(randist);
}

static size_t recombine_rays(intptr_t * ipars, size_t ipars_length, 
		struct r_t ** r, size_t * r_size, size_t preserve){
	size_t r_length = 0;
	intptr_t * index = malloc(ipars_length * sizeof(intptr_t));
	intptr_t j;
	for(j = 0; j < ipars_length; j++) {
		intptr_t ipar = ipars[j];
		if(psys.recomb[ipar] > 0 /*FIXME: use a bigger threshold*/) {
			index[j] = r_length + preserve;
			r_length++;
		}
	}

	if(r_length + preserve > *r_size) {
		size_t old_size = *r_size;
		while(r_length + preserve > *r_size) {
			if(*r_size == 0) *r_size = 1024;
			*r_size *=2;
		}
		*r = realloc(*r, sizeof(struct r_t) * *r_size);
		
		intptr_t i;
		for(i = old_size; i < *r_size; i++) {
			(*r)[i].x = NULL;
			(*r)[i].ipars = NULL;
			(*r)[i].ipars_length = 0;
			(*r)[i].ipars_size = 0;
		}
	}

	if(r_length == 0) {
		free(index);
		return 0;
	}

	for(j = 0; j < ipars_length; j++) {
		intptr_t ipar = ipars[j];
		intptr_t i = index[j];
		if(psys.recomb[ipar] > 0 /*FIXME: use a bigger threshold*/) {
			(*r)[i].s[0] = psys.pos[ipar][0];
			(*r)[i].s[1] = psys.pos[ipar][1];
			(*r)[i].s[2] = psys.pos[ipar][2];
			double dx,dy,dz;
			#pragma omp atmoic
			gsl_ran_dir_3d(rng, &dx, &dy, &dz);
			(*r)[i].dir[0] = dx;
			(*r)[i].dir[1] = dy;
			(*r)[i].dir[2] = dz;
			(*r)[i].Nph = psys.recomb[ipar];
			(*r)[i].freq = 1.0;
			(*r)[i].length = psys.boxsize;
			psys.recomb[ipar] = 0;
		}
	}
	free(index);
	return r_length;
}

static void trace(struct r_t * r, size_t Nr) {
	intptr_t i;
	for(i = 0; i < Nr; i++) {
		size_t ipars_size_new = r[i].ipars_size;
		r[i].ipars_length = rt_trace(r[i].s, r[i].dir, r[i].length, &r[i].ipars, &ipars_size_new);
		if(ipars_size_new > r[i].ipars_size) {
			free(r[i].x);
			r[i].x = malloc(sizeof(struct x_t) * ipars_size_new);
		}
		r[i].ipars_size = ipars_size_new;
	}

}
static size_t deposit(struct r_t * r, size_t Nr, intptr_t **ipars, size_t * ipars_size) {

	intptr_t i;
	for(i = 0; i < Nr; i++) {
		intptr_t j; /*index of the partilce in pars */
		/* clear the deposit of the relavant particles */
		/* this shall loop over all rays */
		for(j = 0; j < r[i].ipars_length; j++) {
			intptr_t ipar = r[i].ipars[j];
			psys.deposit[ipar] = 0.0;
		}
	}

	for(i = 0; i < Nr; i++) {
		intptr_t j; /*index of the partilce in pars */
		/* sort particles by distance from source */
		for(j = 0; j < r[i].ipars_length; j++) {
			intptr_t ipar = r[i].ipars[j];
			r[i].x[j].ipar = ipar;
			float dist = 0.0;
			float proj = 0.0;
			int d;
			for(d = 0; d < 3; d++) {
				float dd = psys.pos[ipar][d] - r[i].s[d];
				dist += dd * dd;
				proj += dd * r[i].dir[d];
			}
			r[i].x[j].d = sqrt(dist);
			r[i].x[j].b = sqrt(fabs(dist - proj * proj));
		}
		qsort(r[i].x, r[i].ipars_length, sizeof(struct x_t), x_t_compare);
	}

	for(i = 0; i < Nr; i++) {
		double Tau = 0.0;
		double TM = r[i].Nph; /*transmission*/
		intptr_t j;
		for(j = 0; j < r[i].ipars_length; j++) {
			intptr_t ipar = r[i].x[j].ipar;
			float b = r[i].x[j].b;
			double sigma = ar_verner(r[i].freq);
			float sml = psys.sml[ipar];
			double NHI = psys.xHI[ipar] * psys.mass[ipar] * C_HMF / U_MPROTON;
			double Ncd = sph_depth(b / sml) / (sml * sml) * NHI;
			float tau = sigma * Ncd;

			double absorb = 0.0;
			if(tau < 0.001) absorb = TM * tau;
			else absorb = TM * (1. - exp(-tau));

			/* make this atomic */
			psys.deposit[ipar] += absorb;
			bitmask_set(active, ipar);

			TM -= absorb;

			if(Tau > 30) {
				r[i].ipars_length = j;
				break;
			}
		}
	}

	size_t ipars_length_est = 0;
	for(i = 0; i < Nr; i++) {
		ipars_length_est += r[i].ipars_length;
	}

	if(ipars_length_est > *ipars_size) {
		*ipars_size = ipars_length_est;
		free(*ipars);
		*ipars = malloc(sizeof(intptr_t) * *ipars_size);
	}

	size_t ipars_length = 0;
	for(i = 0; i < Nr; i++) {
		intptr_t j;
		for(j = 0; j < r[i].ipars_length; j++) {
			intptr_t ipar = r[i].ipars[j];
			if(bitmask_test_and_clear(active, ipar)) {
				(*ipars)[ipars_length] = ipar;
				ipars_length++;
			}
		}
	}
	return ipars_length;
}

static void solve(intptr_t * ipars, size_t ipars_length) {
	Solver * s = solver_new();
	intptr_t j;
	for(j = 0; j < ipars_length; j++) {
		intptr_t ipar = ipars[j];
		solver_evolve(s, ipar);
		psys.deposit[ipar] = 0;
		psys.lasthit[ipar] = psys.tick;
	}
	solver_delete(s);
}

	
static int x_t_compare(const void * p1, const void * p2) {
	const struct x_t * x1 = p1;
	const struct x_t * x2 = p2;
	if(x1->d > x2->d) return 1;
	if(x1->d < x2->d) return -1;
	if(x1->d == x2->d) return 0;
}
static float dist(const float p1[3], const float p2[3]) {
	int d;
	double result = 0.0;
	for (d = 0; d < 3; d++) {
		float dd = fabs(p1[d] - p2[d]);
		result += dd *dd;
	}
	return sqrt(result);
}
