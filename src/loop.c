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

struct r_t {
	float s[3];
	float dir[3];
	double Nph;
	double freq;
	double length;

	int type;
	union {
		intptr_t isrc;
		intptr_t ipar;
	};

	ARRAY_DEFINE_S(x, Xtype)
};

extern PSystem psys;
typedef struct _Solver Solver;
extern Solver * solver_new();
extern int solver_evolve(Solver * s, intptr_t ipar);
extern void solver_delete(Solver * s);
extern float sph_depth(const float r_h);

size_t rt_trace(const float s[3], const float dir[3], const float dist, Xtype ** pars, size_t * size);

static int x_t_compare(const void * p1, const void * p2);
static int r_t_compare(const void * p1, const void * p2);

static void trace();
static void emit_rays();
static void merge_ipars();
static void deposit();
static void update_pars();

ARRAY_DEFINE(r, struct r_t);
ARRAY_DEFINE(ipars, intptr_t);

char * active = NULL;
static struct {
	struct {
		double recomb;
		double src;
		double lost;
	} tick_photon;
	struct {
		double recomb;
		double src;
		double lost;
	} total_photon;
	struct {
		size_t recomb;
		size_t src;
	} total_ray;
	struct {
		size_t recomb;
		size_t src;
	} tick_ray;

	size_t destruct;
	size_t errors;
	size_t count;
	double recomb_pool_photon;
	size_t tick;
} stat;

void init() {
}

void run_epoch() {

	active = bitmask_alloc(psys.npar);
	bitmask_clear_all(active);

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
			stat.total_ray.src += stat.tick_ray.src;
			stat.total_ray.recomb += stat.tick_ray.recomb;
			stat.total_photon.src += stat.tick_photon.src;
			stat.total_photon.lost += stat.tick_photon.lost;
			stat.total_photon.recomb += stat.tick_photon.recomb;

			MESSAGE("-----tick: %lu ---[SRC REC]----", psys.tick);
			MESSAGE("RY %lu/%lu %lu/%lu", 
			        stat.tick_ray.src / stat.tick, stat.total_ray.src, 
			        stat.tick_ray.recomb / stat.tick, stat.total_ray.recomb);
			MESSAGE("PH %le/%le %le/%le %le/%le", 
			        stat.tick_photon.src / stat.tick, stat.total_photon.src,
			        stat.tick_photon.recomb / stat.tick, stat.total_photon.recomb,
			        stat.tick_photon.lost / stat.tick, stat.total_photon.lost);
			MESSAGE("RP: %le", stat.recomb_pool_photon);
			MESSAGE("EV: %lu/%lu", stat.errors, stat.count);
			psystem_stat("T");
			psystem_stat("xHI");
			psystem_stat("ye");
			psystem_stat("ie");
			psystem_stat("recomb");

			psystem_write_output(istep + 1);
			istep++;
			memset(&stat.tick_ray, 0, sizeof(stat.tick_ray));
			memset(&stat.tick_photon, 0, sizeof(stat.tick_photon));
			stat.tick = 0;
		}
		if(psys.tick == psys.epoch->nticks) break;

		stat.tick++;
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
		ARRAY_FREE(r[i].x);
	}
	ARRAY_FREE(r);

	free(active);
}


static void emit_rays() {
	double weights[psys.nsrcs];
	double total_src = 0.0;
	intptr_t i;
	for(i = 0; i < psys.nsrcs; i++) {
		weights[i] = psys.srcs[i].Ngamma_sec * psys.tick_time / U_SEC;
		total_src += weights[i];
	}

	
	gsl_ran_discrete_t * src_ran = NULL;
	gsl_ran_discrete_t * rec_ran = NULL;
	const double branch = total_src / (total_src + stat.recomb_pool_photon);
	for(i = 0; i < r_length; i++) {
		const double sample = gsl_rng_uniform(RNG);
		if(!CFG_DISABLE_2ND_GEN_PHOTONS 
		&& sample > branch) {
			/* recomb ray */
			if(rec_ran == NULL) rec_ran = gsl_ran_discrete_preproc(psys.npar, psys.recomb);
			const intptr_t ipar = gsl_ran_discrete(RNG, rec_ran);
			r[i].ipar = ipar;
			r[i].type = 1;
		} else {
			if(src_ran == NULL) src_ran = gsl_ran_discrete_preproc(psys.nsrcs, weights);
			const intptr_t isrc = gsl_ran_discrete(RNG, src_ran);
			r[i].isrc = isrc;
			/* src ray */
			r[i].type = 0;
		}
	}

	/* sort the rays by type and ipar/isrc so that
     * so that we know how many times a source is sampled */
	qsort(r, r_length, sizeof(struct r_t), r_t_compare);
    
	/* now make the photons in the rays */
	int identical_count = -10000;
	for(i = 0; i < r_length; i++) {
		int d;
		double dx, dy, dz;
		switch(r[i].type) {
			case 0: /* from a src */
			for(d = 0; d < 3; d++) {
				r[i].s[d] = psys.srcs[r[i].isrc].pos[d];
			}
			r[i].Nph = weights[r[i].isrc];
			r[i].freq = spec_gen_freq(psys.srcs[r[i].isrc].specid);
			break;
			case 1: /* from a recombination, aka particle */
			for(d = 0; d < 3; d++) {
				r[i].s[d] = psys.pos[r[i].ipar][d];
			}
			r[i].Nph = psys.recomb[r[i].ipar]; 
			r[i].freq = 1.0;
			break;
			default: ERROR("never each here");
		}
		gsl_ran_dir_3d(RNG, &dx, &dy, &dz);
		r[i].dir[0] = dx;
		r[i].dir[1] = dy;
		r[i].dir[2] = dz;
		r[i].length = 2.0 * psys.boxsize;
	}
	for(i = 0; i <= r_length; i++) {
		if(i == 0) {
			identical_count = 1;
			continue;
		}
		if(i == r_length || r_t_compare(&r[i-1], &r[i])) {
			const double factor = 1.0 / identical_count;
			while(identical_count > 0) {
				r[i - identical_count].Nph *= factor;
				identical_count --;
			}
			identical_count = 1;
		} else {
			identical_count ++;
		}
	}

	/* rays generated, do some stat and purge the recomb pool in psys*/
	for(i = 0; i < r_length; i++) {
		switch(r[i].type) {
			case 0: /* from a src */
				stat.tick_ray.src++;
				stat.tick_photon.src += r[i].Nph;
			break;
			case 1: /* from a recombination, aka particle */
				stat.tick_ray.recomb++;
				stat.tick_photon.recomb += r[i].Nph;
				stat.recomb_pool_photon -= r[i].Nph;
				psys.recomb[r[i].ipar] -= r[i].Nph;
				if(psys.recomb[r[i].ipar] < 0) {
					psys.recomb[r[i].ipar] = 0;
				}
			break;
			default: ERROR("never each here");
		}
	}

	if(src_ran != NULL) gsl_ran_discrete_free(src_ran);
	if(rec_ran != NULL) gsl_ran_discrete_free(rec_ran);
}

static void trace() {
	intptr_t i;
	#pragma omp parallel for private(i)
	for(i = 0; i < r_length; i++) {
		r[i].x_length = rt_trace(r[i].s, r[i].dir, r[i].length, 
			&r[i].x, &r[i].x_size);

		qsort(r[i].x, r[i].x_length, sizeof(Xtype), x_t_compare);
	}

}

static void deposit(){
	intptr_t i;

	const double U_CM2 = U_CM * U_CM;
	const double NH_fac = C_HMF / U_MPROTON;

	for(i = 0; i < r_length; i++) {
		double Tau = 0.0;
		double TM = r[i].Nph; /*transmission*/
		intptr_t j;
		for(j = 0; j < r[i].x_length; j++) {
			const intptr_t ipar = r[i].x[j].ipar;
			const float b = r[i].x[j].b;
			const double sigma = ar_verner(r[i].freq) * U_CM2;
			const float sml = psys.sml[ipar];
			const float sml_inv = 1.0 / sml;
			const double NH = psys.mass[ipar] * NH_fac;
			const double NHI = psys.xHI[ipar] * NH;
			const double Ncd = sph_depth(b * sml_inv) * (sml_inv * sml_inv) * NHI;
			const double tau = sigma * Ncd;

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
			
			double newxHI = psys.xHI[ipar] - absorb / NH;
			if(newxHI < 0.0) newxHI = 0.0;
			if(newxHI > 1.0) newxHI = 1.0;
			
			const double dxHI = newxHI - psys.xHI[ipar];
			const double newye = psys.ye[ipar] - dxHI;

			if(CFG_ISOTHERMAL) {
				const double T = ieye2T(psys.ie[ipar], psys.ye[ipar]);
				psys.ie[ipar] = Tye2ie(T, newye);
			}
			psys.xHI[ipar] = newxHI;
			psys.ye[ipar] = newye;
			
			Tau += tau;
			/* cut off at around 10^-10 */
			if(Tau > 30.0) {
				ARRAY_RESIZE(r[i].x, Xtype, j);
				break;
			}
		}
		stat.tick_photon.lost += TM;
	}
}

static void merge_ipars() {
	intptr_t i;
	/* now merge the list */
	size_t ipars_length_est = 0;
	for(i = 0; i < r_length; i++) {
		ipars_length_est += r[i].x_length;
	}

	ARRAY_ENSURE(ipars, intptr_t, ipars_length_est);
	ARRAY_CLEAR(ipars);

	for(i = 0; i < r_length; i++) {
		intptr_t j;
		for(j = 0; j < r[i].x_length; j++) {
			bitmask_set(active, r[i].x[j].ipar);
		}
	}

	for(i = 0; i < r_length; i++) {
		intptr_t j;
		for(j = 0; j < r[i].x_length; j++) {
			intptr_t ipar = r[i].x[j].ipar;
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
	const double NH_fac = C_HMF / U_MPROTON;
	const double nH_fac = C_HMF / (U_MPROTON / (U_CM * U_CM * U_CM));

	#pragma omp parallel for reduction(+: d1, d2, increase_recomb) private(j) schedule(static)
	for(j = 0; j < ipars_length; j++) {
		const intptr_t ipar = ipars[j];
		Step step = {0};
/*
		if(psys.tick == psys.lasthit[ipar]) {
			ERROR("particle solved twice at one tick");
		}
*/
		const double NH = NH_fac * psys.mass[ipar];
		/* everything multiplied by nH, saving some calculations */
		step.xHI = psys.xHI[ipar];
		step.ye = psys.ye[ipar];
		step.y = psys.ye[ipar] - (1.0 - psys.xHI[ipar]);
		step.nH = nH_fac * psys.rho[ipar];
		step.ie = psys.ie[ipar];
		step.T = ieye2T(psys.ie[ipar], psys.ye[ipar]);

		const double time = (psys.tick - psys.lasthit[ipar]) * psys.tick_time;
		if(!step_evolve(&step, time)) {
			d1++;
		} else {
			const double recphotons = step.dyGH * NH;
			psys.recomb[ipar] += recphotons;
			if(psys.recomb[ipar] < 0) psys.recomb[ipar] = 0;
			increase_recomb += recphotons;
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
	stat.recomb_pool_photon += increase_recomb;
}

static int x_t_compare(const void * p1, const void * p2) {
	const double d1 = ((const Xtype * )p1)->d;
	const double d2 = ((const Xtype * )p2)->d;
	if(d1 < d2) return -1;
	if(d1 > d2) return 1;
	if(d1 == d2) return 0;
}
static int r_t_compare(const void * p1, const void * p2) {
	const int t1 = ((const struct r_t * )p1)->type;
	const int t2 = ((const struct r_t * )p2)->type;
	if(t1 != t2) {
		return t1 - t2;
	}
	const intptr_t i1 = ((const struct r_t * )p1)->ipar;
	const intptr_t i2 = ((const struct r_t * )p2)->ipar;
	/* note that intptr_t is longer than int, truncating may make errors*/
	if(i1 > i2) return -1;
	if(i1 < i2) return 1;
	return 0;
}
