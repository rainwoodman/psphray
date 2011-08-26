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
#ifdef _OPENMP
#include <omp.h>
#endif
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

extern float sph_depth(const float r_h);
extern int step_evolve(Step * step, double time);

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

	struct {
		double src_length;
		double rec_length;
	} ray;

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

	memset(&stat, 0, sizeof(stat));

	intptr_t istep = 0;
	psystem_stat("xHI");
	psystem_stat("ye");
	psystem_stat("yeMET");
	psystem_stat("ie");
	psystem_stat("T");
	psystem_stat("yGrec");
	while(1) {
		if(istep < psys.epoch->output.nsteps && psys.tick == psys.epoch->output.steps[istep]) {
			stat.total_ray.src += stat.tick_ray.src;
			stat.total_ray.recomb += stat.tick_ray.recomb;
			stat.total_photon.src += stat.tick_photon.src;
			stat.total_photon.lost += stat.tick_photon.lost;
			stat.total_photon.recomb += stat.tick_photon.recomb;

			MESSAGE("-----tick: %lu ---[SRC REC]----", psys.tick);
			MESSAGE("RY %lu/%lu(%g) %lu/%lu(%g)", 
			        stat.tick_ray.src / stat.tick, stat.total_ray.src, 
					stat.ray.src_length,
			        stat.tick_ray.recomb / stat.tick, stat.total_ray.recomb,
					stat.ray.rec_length);
			MESSAGE("PH %le/%le %le/%le %le/%le [leak]", 
			        stat.tick_photon.src / stat.tick, stat.total_photon.src,
			        stat.tick_photon.recomb / stat.tick, stat.total_photon.recomb,
			        stat.tick_photon.lost / stat.tick, stat.total_photon.lost);
			MESSAGE("RP: %le, NRAY_MAX %ld", stat.recomb_pool_photon, r_size);
			MESSAGE("EV: %lu/%lu", stat.errors, stat.count);
			psystem_stat("T");
			psystem_stat("xHI");
			psystem_stat("yeMET");
			psystem_stat("ye");
			psystem_stat("ie");
			psystem_stat("yGrec");

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
	intptr_t i;
	for(i = 0; i < psys.nsrcs; i++) {
		weights[i] = psys.srcs[i].Ngamma_sec * (psys.tick - psys.srcs[i].lastemit) * psys.tick_time / U_SEC;
	}

	double max_src_length = -1.0;
	double max_rec_length = -1.0;
	for(i = 0; i < r_length; i++) {
		if(r[i].type == 0 && max_src_length < r[i].length) max_src_length = r[i].length;
		if(r[i].type == 1 && max_rec_length < r[i].length) max_rec_length = r[i].length;
	}
	if(max_src_length < 0 || max_src_length > 2.0 * psys.boxsize) {
		max_src_length = 2.0 * psys.boxsize;
	}
	if(max_rec_length < 0 || max_rec_length > 2.0 * psys.boxsize) {
		max_rec_length = 2.0 * psys.boxsize;
	}

	stat.ray.src_length = max_src_length;
	stat.ray.rec_length = max_rec_length;

	ARRAY_RESIZE(r, struct r_t, psys.epoch->nray);

	gsl_ran_discrete_t * src_ran = gsl_ran_discrete_preproc(psys.nsrcs, weights);
	for(i = 0; i < psys.epoch->nray; i++) {
		const intptr_t isrc = gsl_ran_discrete(RNG, src_ran);
		r[i].isrc = isrc;
		/* src ray */
		r[i].type = 0;
	}

	if(ipars_length > 0) {
		size_t j = psys.epoch->nray;
		double fsum = 0.0;
		double f[ipars_length];
		#pragma omp parallel for private(i) reduction(+: fsum)
		for(i = 0; i < ipars_length; i++) {
			const intptr_t ipar = ipars[i];
			const double rec = psys.yGrec[ipar] * psys_NH(ipar);
			f[i] = psys.yGrec[ipar];
			fsum += f[i];
		}

		gsl_ran_discrete_t * rec_ran = gsl_ran_discrete_preproc(ipars_length, f);

		intptr_t k;
		for(k = 0; k < psys.epoch->nrec; k++) {
			int i = gsl_ran_discrete(RNG, rec_ran);
			ARRAY_RESIZE(r, struct r_t, j + 1);
			r[j].type = 1;
			r[j].ipar = ipars[i];
			j++;
		}
		gsl_ran_discrete_free(rec_ran);
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
			r[i].length = max_src_length;
			break;
			case 1: /* from a recombination, aka particle */
			{
			const intptr_t ipar = r[i].ipar;
			for(d = 0; d < 3; d++) {
				r[i].s[d] = psys.pos[ipar][d];
			}
			r[i].Nph = psys.yGrec[ipar] * psys_NH(ipar);
			r[i].freq = 1.0;
			r[i].length = max_rec_length;
			break;
			}
			default: ERROR("never each here");
		}
		gsl_ran_dir_3d(RNG, &dx, &dy, &dz);
		r[i].dir[0] = dx;
		r[i].dir[1] = dy;
		r[i].dir[2] = dz;
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
				psys.srcs[r[i].isrc].lastemit = psys.tick;
			break;
			case 1: /* from a recombination, aka particle */
			{
				const intptr_t ipar = r[i].ipar;
				stat.tick_ray.recomb++;
				stat.tick_photon.recomb += r[i].Nph;
				stat.recomb_pool_photon -= r[i].Nph;
				psys.yGrec[ipar] -= r[i].Nph / psys_NH(ipar);
				if(psys.yGrec[ipar] < 0) {
					psys.yGrec[ipar] = 0;
				}
			}
			break;
			default: ERROR("never each here");
		}
	}

	if(src_ran != NULL) gsl_ran_discrete_free(src_ran);
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

	const double scaling_fac2_inv = CFG_COMOVING?pow((psys.epoch->redshift + 1),2):1.0;
	const double scaling_fac = CFG_COMOVING?1/(psys.epoch->redshift + 1):1.0;
	#pragma omp parallel for private(i)
	for(i = 0; i < r_length; i++) {
		double Tau = 0.0;
		double TM = r[i].Nph; /*transmission*/
		intptr_t j;
		for(j = 0; j < r[i].x_length; j++) {
			const intptr_t ipar = r[i].x[j].ipar;
			while(bitmask_test_and_set(active, r[i].x[j].ipar)) {
				continue;
			}
			const float b = r[i].x[j].b;
			const double sigma = xs_get(XS_HI, r[i].freq) * U_CM2;
			const float sml = psys.sml[ipar];
			const float sml_inv = 1.0 / sml;
			const double NH = psys_NH(ipar);
			const double xHI = psys_xHI(ipar);
			const double xHII = psys_xHII(ipar);

			const double NHI = xHI * NH;
			const double Ncd = sph_depth(b * sml_inv) * (sml_inv * sml_inv) * NHI * scaling_fac2_inv;
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
			
			const double delta = absorb / NH;

			psys.yGdep[ipar] += delta;

			bitmask_clear(active, r[i].x[j].ipar);
			/* cut off at around 10^-10 */
			if(TM / r[i].Nph < 1e-10) {
				ARRAY_RESIZE(r[i].x, Xtype, j);
				/* point to the next intersection so that a few lines later we can use j - 1*/
				j++;
				break;
			}
		}
		// enlarge the length by a bit 
		// if we terminated the ray sooner then it is safe to do this
		// if the ray terminated too early we shall use a longer length
		// next time.
		r[i].length = fmax(r[i].x[j - 1].d * 2.0,  r[i].x[j - 1].d + C_SPEED_LIGHT * scaling_fac * psys.tick_time);
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
	const double nH_fac = C_HMF / (U_MPROTON / (U_CM * U_CM * U_CM));
	const double scaling_fac3_inv = CFG_COMOVING?pow((psys.epoch->redshift + 1.0), 3):1.0;
	#pragma omp parallel for reduction(+: d1, d2, increase_recomb) private(j) schedule(static)
	for(j = 0; j < ipars_length; j++) {
		const intptr_t ipar = ipars[j];
		const double delta = psys.yGdep[ipar];
		const double xHI = psys_xHI(ipar);
		const double xHII = psys_xHII(ipar);
		if(CFG_ISOTHERMAL) {
			const double T = psys_T(ipar);
			psys_set_lambdaHI(ipar, xHI - delta, xHII + delta);
			psys.ie[ipar] = Tye2ie(T, psys_ye(ipar));
		} else {
			psys_set_lambdaHI(ipar, xHI - delta, xHII + delta);
		}
		psys.yGdep[ipar] = 0.0;
		Step step = {0};
/*
		if(psys.tick == psys.lasthit[ipar]) {
			ERROR("particle solved twice at one tick");
		}
*/
		const double NH = psys_NH(ipar);
		step.yGdep = psys.yGdep[ipar];
		step.lambdaHI = psys.lambdaHI[ipar];
		step.yeMET = psys.yeMET[ipar];
		step.nH = nH_fac * psys.rho[ipar] * scaling_fac3_inv;
		step.ie = psys.ie[ipar];
		step.T = psys_T(ipar);

		const double time = (psys.tick - psys.lasthit[ipar]) * psys.tick_time;
		if(!step_evolve(&step, time)) {
			WARNING("evolve failed: time,T,lambdaHI,y,nH,ie=%g %g %g %g %g %g %g",
				time/U_MYR, step.T, step.lambdaHI, step.yeMET, step.nH, step.ie);
			abort();
			d1++;
		} else {
			psys.yGrec[ipar] += step.dyGrec;
			if(psys.yGrec[ipar] < 0) psys.yGrec[ipar] = 0;
			increase_recomb += step.dyGrec * NH;

			psys.lambdaHI[ipar] = step.lambdaHI;
			psys.ie[ipar] += step.ie;
			psys.lasthit[ipar] = psys.tick;

			if(CFG_ISOTHERMAL) {
				psys.ie[ipar] = Tye2ie(step.T, psys_ye(ipar));
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
