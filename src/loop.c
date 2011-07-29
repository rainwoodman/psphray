#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <bitmask.h>
#include <array.h>
#include <messages.h>

#include "config.h"
#include "psystem.h"
#include "reader.h"

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
	double length;

	intptr_t ires;

	ARRAY_DEFINE_S(ipars, intptr_t)
	ARRAY_DEFINE_S(x, struct x_t)
};

struct res_t {
	double Nph;
	intptr_t isrc;
	intptr_t ipar;
	int type;
};

extern PSystem psys;
typedef struct _Solver Solver;
extern Solver * solver_new();
extern int solver_evolve(Solver * s, intptr_t ipar);
extern void solver_delete(Solver * s);
extern float sph_depth(const float r_h);

size_t rt_trace(const float s[3], const float dir[3], const float dist, intptr_t ** pars, size_t * size);

static int x_t_compare(const void * p1, const void * p2);
static float dist(const float p1[3], const float p2[3]);
static void trace();
static void make_photons();
static void emit_rays();
static void merge_ipars();
static void deposit();
static void update_pars();

gsl_rng * rng = NULL;;

ARRAY_DEFINE(res, struct res_t);
ARRAY_DEFINE(r, struct r_t);
ARRAY_DEFINE(ipars, intptr_t);

char * active = NULL;
static struct {
	struct {
		double recomb;
		double source;
		double lost;
	} ray_photon;
	struct {
		size_t recomb;
		size_t source;
	} ray_count;
	size_t destruct;
	size_t errors;
	size_t count;
	size_t full_scans;
	double total_recomb;
} stat;

void init() {
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, 123456);
}

void run() {

	active = bitmask_alloc(psys.npar);
	bitmask_clear_all(active);
	ARRAY_ENSURE(res, struct res_t, psys.nsrcs);

	ARRAY_RESIZE(r, struct r_t, psys.epoch->nray);

	memset(&stat, 0, sizeof(stat));

	intptr_t istep = 0;
	psystem_stat("xHI");
	psystem_stat("ye");
	psystem_stat("ie");
	psystem_stat("T");
	psystem_stat("recomb");
	while(1) {
		if(istep < psys.epoch->output.nsteps && psys.tick == psys.epoch->output.steps[istep]) {
			MESSAGE("tick: %lu, full recomb scans = %lu", psys.tick, stat.full_scans);
			MESSAGE("Rays: %lu(rec) %lu(src)", stat.ray_count.recomb, stat.ray_count.source);
			MESSAGE("Ray photons: %le(rec) %le(src) %le(lost)", stat.ray_photon.recomb, stat.ray_photon.source, stat.ray_photon.lost);
			MESSAGE("destructive: %lu, evolve error %lu/%lu", stat.destruct, stat.errors, stat.count);
			MESSAGE("recombine pool %le RES pool %lu", stat.total_recomb, res_size);
			psystem_stat("xHI");
			psystem_stat("ye");
			psystem_stat("ie");
			psystem_stat("recomb");

			psystem_write_output(istep + 1);
			istep++;
		}
		if(psys.tick == psys.epoch->nticks) break;

		make_photons();
		psys.tick++;

		emit_rays();

		/* trace rays */
		trace();

		/* deposit photons */
		deposit();

		merge_ipars();

		update_pars();
	}

	ARRAY_FREE(ipars);
	intptr_t i;
	for(i = 0 ;i < r_length; i++) {
		ARRAY_FREE(r[i].ipars);
		ARRAY_FREE(r[i].x);
	}
	ARRAY_FREE(r);

	ARRAY_FREE(res);
	free(active);
}

static void make_photons(){
	intptr_t i;

	ARRAY_RESIZE(res, struct res_t, psys.nsrcs);

	double total_src = 0.0;
	for(i = 0; i < psys.nsrcs; i++) {
		res[i].Nph += psys.srcs[i].Ngamma_sec * psys.tick_time / U_SEC;
		res[i].type = 0;
		res[i].isrc = i;
		total_src += res[i].Nph;
	}

	if(CFG_DISABLE_2ND_GEN_PHOTONS) return;

	if(stat.total_recomb > 0.1 * total_src) {
		/* if there are too many recomb photons, do a full recombineation scan*/
		stat.full_scans++;
		intptr_t ipar;
		for(ipar = 0; ipar < psys.npar; ipar++) {
			if(psys.recomb[ipar] > 0.0 /*FIXME: use a bigger threshold*/) {
				struct res_t * t = ARRAY_APPEND(res, struct res_t);
				t->Nph = psys.recomb[ipar];
				t->type = 1;
				t->ipar = ipar;
			}
		}
	} else {
		for(i = 0; i < ipars_length; i++) {
			intptr_t ipar = ipars[i];
			if(psys.recomb[ipar] > 0.0 /*FIXME: use a bigger threshold*/) {
				struct res_t * t = ARRAY_APPEND(res, struct res_t);
				t->Nph = psys.recomb[ipar];
				t->type = 1;
				t->ipar = ipar;
			}
		}
	}
}

static void emit_rays() {
	double weights[res_length];
	intptr_t i;
	size_t raycount[res_length];
	for(i = 0; i < res_length; i++) {
		switch(res[i].type) {
			case 0:
			weights[i] = res[i].Nph;
			break;
			case 1:
			weights[i] = res[i].Nph;
		}
		raycount[i] = 0;
	}

	gsl_ran_discrete_t * randist = gsl_ran_discrete_preproc(res_length, weights);

	for(i = 0; i < r_length; i++) {
		intptr_t ires = gsl_ran_discrete(rng, randist);
		r[i].ires = ires;
		raycount[ires]++;
	}
	for(i = 0; i < r_length; i++) {
		intptr_t ires = r[i].ires;
		int d;
		double dx, dy, dz;
		switch(res[ires].type) {
			case 0: /* from a src */
			for(d = 0; d < 3; d++) {
				r[i].s[d] = psys.srcs[res[ires].isrc].pos[d];
			}
			break;
			case 1: /* from a recombination, aka particle */
			for(d = 0; d < 3; d++) {
				r[i].s[d] = psys.pos[res[ires].ipar][d];
			}
			break;
			default: ERROR("never each here");
		}
		gsl_ran_dir_3d(rng, &dx, &dy, &dz);
		r[i].dir[0] = dx;
		r[i].dir[1] = dy;
		r[i].dir[2] = dz;
		r[i].freq = 1.0;
		r[i].Nph = res[ires].Nph / raycount[ires];
		r[i].length = 2.0 * psys.boxsize;
	}
	for(i = 0; i < r_length; i++) {
		intptr_t ires = r[i].ires;
		res[ires].Nph = 0.0;
		switch(res[ires].type) {
			case 0: /* from a src */
				stat.ray_count.source++;
				stat.ray_photon.source += r[i].Nph;
			break;
			case 1: /* from a recombination, aka particle */
				stat.ray_count.recomb++;
				stat.ray_photon.recomb += r[i].Nph;
				stat.total_recomb -= r[i].Nph;
				psys.recomb[res[ires].ipar] -= r[i].Nph;
			break;
			default: ERROR("never each here");
		}
	}
	gsl_ran_discrete_free(randist);
}

static void trace() {
	intptr_t i;
	#pragma omp parallel for private(i)
	for(i = 0; i < r_length; i++) {
		r[i].ipars_length = rt_trace(r[i].s, r[i].dir, r[i].length, 
			&r[i].ipars, &r[i].ipars_size);

		ARRAY_RESIZE(r[i].x, struct x_t, r[i].ipars_length);

		intptr_t j;
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

}

static void deposit(){
	intptr_t i;

	for(i = 0; i < r_length; i++) {
		double Tau = 0.0;
		double TM = r[i].Nph; /*transmission*/
		intptr_t j;
		for(j = 0; j < r[i].ipars_length; j++) {
			intptr_t ipar = r[i].x[j].ipar;
			float b = r[i].x[j].b;
			double sigma = ar_verner(r[i].freq) * U_CM * U_CM;
			float sml = psys.sml[ipar];
			double NH = psys.mass[ipar] * C_HMF / U_MPROTON;
			double NHI = psys.xHI[ipar] * NH;
			double Ncd = sph_depth(b / sml) / (sml * sml) * NHI;
			double tau = sigma * Ncd;

			double absorb = 0.0;
			if(tau < 0.001) absorb = TM * tau;
			else if(tau > 100) absorb = TM;
			else absorb = TM * (1. - exp(-tau));

		//	MESSAGE("b = %g transimt = %g absorb = %g Tau=%g", b, TM,absorb, Tau);
			/* make this atomic */
			if(absorb > NHI) {
				stat.destruct++;
				absorb = NHI;
				TM -= NHI;
			} else {
				TM *= exp(-tau);
			}
			
			double T = ieye2T(psys.ie[ipar], psys.ye[ipar]);
			double newxHI = psys.xHI[ipar] - absorb / NH;
			if(newxHI < 0.0) newxHI = 0.0;
			if(newxHI > 1.0) newxHI = 1.0;
			
			double dxHI = newxHI - psys.xHI[ipar];
			psys.ye[ipar] -= dxHI;
			psys.xHI[ipar] = newxHI;
			if(CFG_ISOTHERMAL) {
				psys.ie[ipar] = Tye2ie(T, psys.ye[ipar]);
			}
			Tau += tau;
			/* cut off at around 10^-10 */
			if(Tau > 30.0) {
				ARRAY_RESIZE(r[i].ipars, struct r_t, j);
				break;
			}
		}
		stat.ray_photon.lost += TM;
	}
}

static void merge_ipars() {
	intptr_t i;
	/* now merge the list */
	size_t ipars_length_est = 0;
	for(i = 0; i < r_length; i++) {
		ipars_length_est += r[i].ipars_length;
	}

	ARRAY_ENSURE(ipars, intptr_t, ipars_length_est);
	ARRAY_CLEAR(ipars);

	for(i = 0; i < r_length; i++) {
		intptr_t j;
		for(j = 0; j < r[i].ipars_length; j++) {
			bitmask_set(active, r[i].ipars[j]);
		}
	}

	for(i = 0; i < r_length; i++) {
		intptr_t j;
		for(j = 0; j < r[i].ipars_length; j++) {
			intptr_t ipar = r[i].ipars[j];
			if(bitmask_test_and_clear(active, ipar)) {
				* (ARRAY_APPEND(ipars, intptr_t)) = ipar;
			}
		}
	}
}

static void update_pars() {
	intptr_t j;
	size_t d1 = 0, d2 = 0;
	double increase_recomb = 0;
	#pragma omp parallel for reduction(+: d1, d2, increase_recomb) private(j) schedule(dynamic, 10)
	for(j = 0; j < ipars_length; j++) {
		intptr_t ipar = ipars[j];
		Step step = {0};
		if(psys.tick == psys.lasthit[ipar]) {
			ERROR("particle solved twice at one tick");
		}
		double NH = C_HMF * psys.mass[ipar] / U_MPROTON;
		/* everything multiplied by nH, saving some calculations */
		step.xHI = psys.xHI[ipar];
		step.ye = psys.ye[ipar];
		step.y = psys.ye[ipar] - (1.0 - psys.xHI[ipar]);
		step.nH = C_HMF * psys.rho[ipar] / (U_MPROTON / (U_CM * U_CM * U_CM));
		step.ie = psys.ie[ipar];
		step.T = ieye2T(psys.ie[ipar], psys.ye[ipar]);

		double time = (psys.tick - psys.lasthit[ipar]) * psys.tick_time;
		if(!step_evolve(&step, time)) {
			psystem_stat("T");
			d1++;
		} else {
			psys.recomb[ipar] += step.dyGH * NH;
			increase_recomb += step.dyGH * NH;
			psys.xHI[ipar] += step.dxHI;
			psys.ye[ipar] += step.dye;
			psys.ie[ipar] += step.die;
			psys.lasthit[ipar] = psys.tick;

			if(CFG_ISOTHERMAL) {
				psys.ie[ipar] = Tye2ie(step.T, psys.ye[ipar]);
			}
		}
		d2++;
	}
	stat.errors += d1;
	stat.count += d2;
	stat.total_recomb += increase_recomb;
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
