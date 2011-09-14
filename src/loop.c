#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#define BITMASK_SHUFFLE
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

static inline FILE * fopen_printf(const char * fmt, char * mode, ...) {
	va_list va;
	va_start(va, mode);
	char * str = NULL;
	vasprintf(&str, fmt, va);
	FILE * rt = fopen(str, mode);
	free(str);
	va_end(va);
	return rt;
}
static void maybe_write_particle(const intptr_t ipar, FILE * fp, const char * fmt, ...) {
	if(CFG_DUMP_HOTSPOTS) {
		va_list va;
		va_start(va, fmt);
		intptr_t i;
		if(psys.flag[ipar] & PF_HOTSPOT) {
				vfprintf(fp, fmt, va);
		}
		va_end(va);
	}
}
extern PSystem psys;

extern double sph_depth(const double r_h);
extern double sph_Wh3(const double r_h);

extern int step_evolve(Step * step);

size_t rt_trace(const float s[3], const float dir[3], const float dist, Xtype ** pars, size_t * size);

static int x_t_compare(const void * p1, const void * p2);
static int r_t_compare(const void * p1, const void * p2);

static void trace();
static void emit_rays();
static void merge_pars();
static void deposit();
static void update_pars();

ARRAY_DEFINE(r, struct r_t);
ARRAY_DEFINE(x, Xtype);

static double MAX_SRC_RAY_LENGTH;
static double MAX_REC_RAY_LENGTH;

bitmask_t * active = NULL;
static struct {
	struct {
		size_t total;
		size_t subtotal;
	} src_ray_count;
	struct {
		size_t total;
		size_t subtotal;
	} rec_ray_count;
	struct {
		double total;
		double subtotal;
	} src_photon_count;
	struct {
		double total;
		double subtotal;
	} rec_photon_count;

/* total number of photon travel out of box/ remain in pretruncated ray*/
	double lost_photon_count_sum;
/* total number of recombination photons */
	double rec_photon_count_sum;

	size_t saturated_deposit_count;
	size_t gsl_error_count;
	size_t evolve_count;
	size_t tick_subtotal;

	int * hits;
	FILE * parlogfile;
	FILE * hitlogfile;
} stat;

void init() {
}

void run_epoch() {

	active = bitmask_alloc(psys.npar);
	bitmask_clear_all(active);

	memset(&stat, 0, sizeof(stat));

	stat.hits = calloc(psys.npar, sizeof(int));

	intptr_t istep = 0;

	if(CFG_DUMP_HOTSPOTS) {
		stat.parlogfile = fopen_printf("parlogfile-%03d", "w", psys.epoch - EPOCHS);
		stat.hitlogfile = fopen_printf("hitlogfile-%03d", "w", psys.epoch - EPOCHS);
	}

	psystem_stat("xHI");
	psystem_stat("ye");
	psystem_stat("heat");
	psystem_stat("ie");
	psystem_stat("T");
	psystem_stat("yGrec");
	while(1) {
		if(istep < psys.epoch->output.nsteps && psys.tick == psys.epoch->output.steps[istep]) {
			stat.src_ray_count.total += stat.src_ray_count.subtotal;
			stat.rec_ray_count.total += stat.rec_ray_count.subtotal;
			stat.src_photon_count.total += stat.src_photon_count.subtotal;
			stat.rec_photon_count.total += stat.rec_photon_count.subtotal;

			MESSAGE("-----tick: %lu -------", psys.tick);
			MESSAGE("SRC(subtot): Ray %lu Photon %g Length %g", 
				stat.src_photon_count.subtotal, stat.src_ray_count.subtotal, MAX_SRC_RAY_LENGTH);
			MESSAGE("REC(subtot): Ray %lu Photon %g Length %g", 
				stat.rec_photon_count.subtotal, stat.rec_ray_count.subtotal, MAX_REC_RAY_LENGTH);
			MESSAGE("PH EMISSION: Src %g Rec %g",
				stat.src_photon_count.total, stat.rec_photon_count.total);
			MESSAGE("PH STORAGE : Lost %g Rec %g ", 
				stat.lost_photon_count_sum, stat.rec_photon_count_sum);
			MESSAGE("Evolve     : Error %lu Total %lu", 
				stat.gsl_error_count, stat.evolve_count);

			stat.src_ray_count.subtotal = 0;
			stat.rec_ray_count.subtotal = 0;
			stat.src_photon_count.subtotal = 0;
			stat.rec_photon_count.subtotal = 0;

			psystem_stat("T");
			psystem_stat("xHI");
			psystem_stat("heat");
			psystem_stat("ye");
			psystem_stat("ie");
			psystem_stat("yGrec");

			psystem_write_output(istep + 1);
			FILE * fp = fopen_printf("hits-%03d", "w", istep);
			fwrite(stat.hits, sizeof(int), psys.npar, fp);
			memset(stat.hits, 0, sizeof(int) * psys.npar);
			fclose(fp);
			istep++;
			stat.tick_subtotal = 0;
		}
		if(psys.tick == psys.epoch->nticks) break;

		stat.tick_subtotal++;
		psys.tick++;

		emit_rays();

		/* trace rays */
		trace();

		/* deposit photons */
		deposit();

		merge_pars();

		update_pars();
	}

	ARRAY_FREE(x);
	intptr_t i;
	for(i = 0 ;i < r_length; i++) {
		ARRAY_FREE(r[i].x);
	}
	ARRAY_FREE(r);

	free(active);
	free(stat.hits);
	if(CFG_DUMP_HOTSPOTS) {
		fclose(stat.parlogfile);
		fclose(stat.hitlogfile);
	}
}


static void emit_rays() {
	double weights[psys.nsrcs];
	intptr_t i;
	for(i = 0; i < psys.nsrcs; i++) {
		/* treat two types the same essentially because they are both Ngamma_sec*/
		if(psys.srcs[i].type == PSYS_SRC_POINT) {
			weights[i] = psys.srcs[i].Ngamma_sec * (psys.tick - psys.srcs[i].lastemit) * psys.tick_time / U_SEC;
		} else if(psys.srcs[i].type == PSYS_SRC_PLANE) {
			weights[i] = psys.srcs[i].Ngamma_sec * (psys.tick - psys.srcs[i].lastemit) * psys.tick_time / U_SEC;
		}
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

	MAX_SRC_RAY_LENGTH = max_src_length;
	MAX_REC_RAY_LENGTH = max_rec_length;

	ARRAY_RESIZE(r, struct r_t, psys.epoch->nray);

	gsl_ran_discrete_t * src_ran = gsl_ran_discrete_preproc(psys.nsrcs, weights);
	for(i = 0; i < psys.epoch->nray; i++) {
		const intptr_t isrc = gsl_ran_discrete(RNG, src_ran);
		r[i].isrc = isrc;
		/* src ray */
		r[i].type = 0;
	}

	if(!CFG_ON_THE_SPOT && psys.epoch->nrec && x_length > 0) {
		size_t j = psys.epoch->nray;
		double fsum = 0.0;
		double f[x_length];
		#pragma omp parallel for private(i) reduction(+: fsum)
		for(i = 0; i < x_length; i++) {
			const intptr_t ipar = x[i].ipar;
			const double rec = psys.yGrec[ipar] * psys_NH(ipar);
			f[i] = psys.yGrec[ipar];
			fsum += f[i];
		}

		gsl_ran_discrete_t * rec_ran = gsl_ran_discrete_preproc(x_length, f);

		intptr_t k;
		for(k = 0; k < psys.epoch->nrec; k++) {
			int i = gsl_ran_discrete(RNG, rec_ran);
			ARRAY_RESIZE(r, struct r_t, j + 1);
			r[j].type = 1;
			r[j].ipar = x[i].ipar;
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
			{
			const intptr_t isrc = r[i].isrc;
			if(psys.srcs[isrc].type == PSYS_SRC_POINT) {
				for(d = 0; d < 3; d++) {
					r[i].s[d] = psys.srcs[isrc].pos[d];
				}
				gsl_ran_dir_3d(RNG, &dx, &dy, &dz);
				r[i].dir[0] = dx;
				r[i].dir[1] = dy;
				r[i].dir[2] = dz;
			} else if(psys.srcs[isrc].type == PSYS_SRC_PLANE) {
				const double u = gsl_ran_flat(RNG, -1, 1);
				const double v = gsl_ran_flat(RNG, -1, 1);
				for(d = 0; d < 3; d++) {
					r[i].dir[d] = psys.srcs[isrc].dir[d];
					r[i].s[d] = psys.srcs[isrc].pos[d] + 
						u * psys.srcs[isrc].a[d] * psys.srcs[isrc].radius + 
						v * psys.srcs[isrc].b[d] * psys.srcs[isrc].radius;
				}
				
			}
			r[i].Nph = weights[isrc];
			r[i].freq = spec_gen_freq(psys.srcs[isrc].specid);
			r[i].length = max_src_length;
			}
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
			gsl_ran_dir_3d(RNG, &dx, &dy, &dz);
			r[i].dir[0] = dx;
			r[i].dir[1] = dy;
			r[i].dir[2] = dz;
			break;
			}
			default: ERROR("never each here");
		}
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
				stat.src_ray_count.subtotal++;
				stat.src_photon_count.subtotal += r[i].Nph;
				psys.srcs[r[i].isrc].lastemit = psys.tick;
			break;
			case 1: /* from a recombination, aka particle */
			{
				const intptr_t ipar = r[i].ipar;
				stat.rec_ray_count.subtotal++;
				stat.rec_photon_count.subtotal += r[i].Nph;

				stat.rec_photon_count_sum -= r[i].Nph;
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

	const double scaling_fac = CFG_COMOVING?1/(psys.epoch->redshift + 1):1.0;
	const double scaling_fac2_inv = 1.0 / (scaling_fac * scaling_fac);
	#pragma omp parallel for private(i)
	for(i = 0; i < r_length; i++) {
		double Tau = 0.0;
		double TM = r[i].Nph; /*transmission*/
		intptr_t j;
		const double sigma = xs_get(XS_HI, r[i].freq) * U_CM2;
		for(j = 0; j < r[i].x_length; j++) {
			const intptr_t ipar = r[i].x[j].ipar;
			const double b = r[i].x[j].b;
			const double sml = psys.sml[ipar];
			const double sml_inv = 1.0 / sml;

			int c = 0;
			while(bitmask_test_and_set(active, ipar)) {
				c++;
				if( c > 1000000) {
					WARNING("Dead lock on particle %ld", ipar);
					c = 0;
					continue;
				}
				continue;
			}

			const double NH = psys_NH(ipar);
			const double xHI = psys_xHI(ipar);

			const double NHI = (xHI - psys.yGdep[ipar]) * NH;
			const double Ncd = sph_depth(b * sml_inv) * (sml_inv * sml_inv) * NHI * scaling_fac2_inv;
			double tau = sigma * Ncd;
			if(tau < 0.0) tau = 0.0;
		//	MESSAGE("b = %g transimt = %g absorb = %g Tau=%g", b, TM,absorb, Tau);

			double absorb = - TM * (expm1(-tau));
			if(absorb > NHI) {
				stat.saturated_deposit_count ++;
				absorb = NHI;
				TM -= absorb;
			} else {
				TM *= exp(-tau);

			}
			
			const double delta = absorb / NH;

			psys.yGdep[ipar] += delta;
			bitmask_clear(active, ipar);

			#pragma omp atomic
			psys.heat[ipar] += C_H_PER_MASS * delta * (r[i].freq - 1) * U_RY_ENG;

			#pragma omp atomic
			stat.hits[ipar]++;

			maybe_write_particle(ipar, stat.hitlogfile, 
						"%lu %ld %g %g %g %g %g %g %g %g %g\n",
						psys.tick, ipar, xHI, b, sml, Ncd, sph_depth(b / sml), sigma, NHI, absorb, delta);

			/* cut off at around 10^-10 */
			if(TM / r[i].Nph < 1e-10) {
				/* point to the next intersection so that a few lines later we can use j - 1*/
				j++;
				break;
			}
		}
		// enlarge the length by a bit 
		// if we terminated the ray sooner then it is safe to do this
		// if the ray terminated too early we shall use a longer length
		// next time.
		ARRAY_RESIZE(r[i].x, Xtype, j);
		r[i].length = fmax(r[i].x[j - 1].d * 2.0,  r[i].x[j - 1].d + C_SPEED_LIGHT / scaling_fac * psys.tick_time);
		stat.lost_photon_count_sum += TM;
	}
}

static void merge_pars() {
	intptr_t i;
	/* now merge the list */
	size_t x_length_est = 0;
	for(i = 0; i < r_length; i++) {
		x_length_est += r[i].x_length;
	}

	ARRAY_ENSURE(x, Xtype, x_length_est);
	ARRAY_CLEAR(x);

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
				* (ARRAY_APPEND(x, Xtype)) = r[i].x[j];
			}
		}
	}
}

static void update_pars() {
	intptr_t j;
	size_t d1 = 0, d2 = 0;
	double increase_recomb = 0;
	const double nH_fac = C_HMF / (U_MPROTON / (U_CM * U_CM * U_CM));
	const double scaling_fac = CFG_COMOVING?1/(psys.epoch->redshift + 1.0):1.0;
	const double scaling_fac3_inv = 1.0/(scaling_fac * scaling_fac * scaling_fac);
	#pragma omp parallel for reduction(+: d1, d2, increase_recomb) private(j) schedule(static)
	for(j = 0; j < x_length; j++) {
		const intptr_t ipar = x[j].ipar;
		const double delta = psys.yGdep[ipar];
		const double xHI = psys_xHI(ipar);
		const double xHII = psys_xHII(ipar);
		const double sml_inv = 1.0 / psys.sml[ipar];
		const double rho = psys.rho[ipar];
		
		Step step = {0};

		const double NH = psys_NH(ipar);
#if 1
		step.yGdep = psys.yGdep[ipar];
		step.lambdaHI = psys.lambdaHI[ipar];
#else
		step.lambdaHI = lambdaHI_from_xHI_xHII(xHI - psys.yGdep[ipar], xHII + psys.yGdep[ipar]);
		step.yGdep = 0;
#endif
		step.yeMET = psys.yeMET[ipar];
		step.nH = nH_fac * rho * scaling_fac3_inv;
		step.ie = psys.ie[ipar];
		step.T = psys_T(ipar);
		step.heat = psys.heat[ipar];

		step.time = (psys.tick - psys.lasthit[ipar]) * psys.tick_time;
		maybe_write_particle(ipar, stat.parlogfile, "%lu %lu %g %g %g %g %g %g %g\n",
				psys.tick, ipar, step.time, step.T, step.lambdaHI, step.nH, step.yGdep, step.ie, step.heat);

		if(!step_evolve(&step)) {
			double xHI, xHII;
			lambdaHI_to_xHI_xHII(step.lambdaHI, xHI, xHII);
			WARNING("evolve failed: "
				"time = %g\n"
				"step.time = %g; \n"
				"step.T = %g; \n"
				"step.lambdaHI=%g; \n"
				"step.yeMET=%g;\n"
				"step.nH=%g;\n"
				"step.yGdep=%g;\n"
				"step.ie=%g;\n"
				"step.heat=%g;\n"
				"xHI=%g\n"
				"xHII=%g\n",
				step.time/U_MYR, step.time, step.T, 
				step.lambdaHI, 
				step.yeMET, step.nH, step.yGdep, step.ie, step.heat,
				xHI, xHII
				);
			d1++;
		} else {
			psys.yGrec[ipar] += step.dyGrec;
			if(psys.yGrec[ipar] < 0) psys.yGrec[ipar] = 0;
			increase_recomb += step.dyGrec * NH;

			psys.lambdaHI[ipar] = step.lambdaHI;
			if(!CFG_ADIABATIC)
				psys.ie[ipar] = step.ie;

			psys.lasthit[ipar] = psys.tick;

			if(CFG_ISOTHERMAL) {
				psys.ie[ipar] = Tye2ie(step.T, psys_ye(ipar));
			}
			psys.yGdep[ipar] = 0.0;
			psys.heat[ipar] = 0.0;
		}
		d2++;
	}
	stat.gsl_error_count += d1;
	stat.evolve_count += d2;
	stat.rec_photon_count_sum += increase_recomb;
}

static int x_t_compare(const void * p1, const void * p2) {
	const double d1 = ((const Xtype * )p1)->d;
	const double d2 = ((const Xtype * )p2)->d;
    return (d1 > d2) - (d2 > d1);
}
static int r_t_compare(const void * p1, const void * p2) {
	const int t1 = ((const struct r_t * )p1)->type;
	const int t2 = ((const struct r_t * )p2)->type;
	const intptr_t i1 = ((const struct r_t * )p1)->ipar;
	const intptr_t i2 = ((const struct r_t * )p2)->ipar;
	/* note that intptr_t is longer than int, truncating may make errors*/
	return (((signed int)((t1 > t2) - (t2 > t1))) * 8 ) + ((i1 > i2) - (i2 > i1));
}
