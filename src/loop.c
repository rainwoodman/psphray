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
#include "stat.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <omp.h>
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

extern int step_evolve(void * control, Step * step);
extern void * step_evolve_prepare();
extern void step_evolve_free(void * control);

size_t rt_trace(const float s[3], const float dir[3], const float dist, Xtype ** pars, size_t * size);

static int x_t_compare(const void * p1, const void * p2);
static int r_t_compare(const void * p1, const void * p2);
static int intptr_compare(const void * p1, const void * p2);

static void trace();
static void emit_rays();
static void merge_pars();
static void deposit();
static void update_pars();

ARRAY_DEFINE(r, struct r_t);
ARRAY_DEFINE(active_set, intptr_t);
static size_t x_src_length;

static double MAX_SRC_RAY_LENGTH = 0.0;
static double MAX_REC_RAY_LENGTH = 0.0;

bitmask_t * active = NULL;

void init() {
}

void run_epoch() {

	active = bitmask_alloc(psys.npar);
	bitmask_clear_all(active);
	memset(&stat, 0, sizeof(stat));

	double t0 = omp_get_wtime();
	intptr_t istep = 0;

	stat_restart();

	psystem_stat("xHI");
	psystem_stat("ye");
	psystem_stat("heat");
	psystem_stat("ie");
	psystem_stat("T");
	while(1) {
		if(istep < psys.epoch->output.nsteps && psys.tick == psys.epoch->output.steps[istep]) {
			MESSAGE("-----tick: %lu -------", psys.tick);
			stat_subtotal();

			psystem_stat("T");
			psystem_stat("xHI");
			psystem_stat("xHeI");
			psystem_stat("heat");
			psystem_stat("ye");
			psystem_stat("ie");
			psystem_stat("yGdepHI");
			psystem_stat("yGdepHeI");
			psystem_stat("yGdepHeII");
			psystem_stat("yGrecHII");
			psystem_stat("yGrecHeII");
			psystem_stat("yGrecHeIII");

			psystem_write_output(istep + 1);

			istep++;
			MESSAGE("-----end of tick: %lu -------", psys.tick);
		}
		if(psys.tick == psys.epoch->nticks) break;

		stat.tick_subtotal++;
		psys.tick++;

		double t0 = omp_get_wtime();
		emit_rays();

		/* trace rays */
		trace();

		/* deposit photons */
		if(!CFG_TRACE_ONLY)
		deposit(); 

		merge_pars();

		if(!CFG_TRACE_ONLY)
		update_pars(); 
		stat.total_time += omp_get_wtime() - t0;
	}

	ARRAY_FREE(active_set);
	intptr_t i;
	for(i = 0 ;i < r_length; i++) {
		ARRAY_FREE(r[i].x);
	}
	ARRAY_FREE(r);

	free(active);
	stat_stop();
}


static void emit_rays() {
	const double t0 = omp_get_wtime();
	const double scaling_fac = CFG_COMOVING?1/(psys.epoch->redshift + 1):1.0;
	double weights[psys.nsrcs];
	intptr_t i;
	double max_rec_length = psys.boxsize * 2 * 0.1;

	MAX_REC_RAY_LENGTH = max_rec_length;

	psystem_weight_srcs(weights);

	gsl_ran_discrete_t * src_ran = gsl_ran_discrete_preproc(psys.nsrcs, weights);

	double weightsum = 0.0;
	for(i = 0; i < psys.nsrcs; i++) {
		weightsum += weights[i];
	}
	size_t nray = weightsum / psys.epoch->packet_size;
	size_t nrec = (CFG_ON_THE_SPOT?0:(psys.epoch->nrec * (CFG_H_ONLY?1:3)));
	if(nrec < 0) nrec = 1; /* just to make sure the if 10 lines below works properly */
	if(nray > psys.epoch->nray * 2) {
		ERROR("needs more rays %ld < %ld, psys.epoch->packet_size = %g", psys.epoch->nray, nray, psys.epoch->packet_size);
	}
	ARRAY_RESIZE(r, struct r_t, psys.epoch->nray + nrec);

	r_length = 0;
	for(i = 0; i < psys.epoch->nray; i++) {
		const intptr_t isrc = gsl_ran_discrete(RNG, src_ran);
		r[r_length].isrc = isrc;
		/* src ray */
		r[r_length].type = -1;
		r_length ++;
	}

	if(nrec && psys.tick > 1) {
		int species;
		for(species = 0; species < (CFG_H_ONLY?1:3); species++) {
			if(psys.epoch->nrec > 0) {
				double * weights = malloc(sizeof(double) * x_src_length);
				intptr_t j;
				for(j = 0; j < x_src_length; j++) {
					intptr_t ipar = active_set[j];
					double x = psys_xHI(ipar);
					if(x == 0) x = 1e-10;
					weights[j] = psys.yGrec[species][ipar] / x;
				}
				gsl_ran_discrete_t * ran_rec = gsl_ran_discrete_preproc(x_src_length, weights);
				free(weights);
				intptr_t k;
				for(k = 0; k < psys.epoch->nrec; k++) {
					intptr_t ix = gsl_ran_discrete(RNG, ran_rec);
					struct r_t * p = ARRAY_APPEND_REUSE(r, struct r_t);
					p->type = species;
					p->ipar = active_set[ix];
				}

				gsl_ran_discrete_free(ran_rec);
			} else {
				intptr_t j;
				for(j = 0; j < x_src_length; j++) {
					intptr_t ipar = active_set[j];
					if(psys.yGrec[species][ipar] > CFG_RECOMBINE_THRESHOLD) {
						struct r_t * p = ARRAY_APPEND_REUSE(r, struct r_t);
						p->type = species;
						p->ipar = ipar;
					}
				}
			}
		}
	}

	/* sort the rays by type and ipar/isrc so that
     * so that we know how many times a source is sampled */
	qsort(r, r_length, sizeof(struct r_t), r_t_compare);
    
	/* now make the photons in the rays */
	for(i = 0; i < r_length; i++) {
		int d;
		double dx, dy, dz;
		switch(r[i].type) {
			case -1: /* from a src */
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
			r[i].length = psys.boxsize * 2; //psys.srcs[isrc].ray_length_hint;
			if(r[i].freq > 100000 * U_EV) {
				ERROR("freq strangely large");
			}
			}
			break;
			case 0: /* from a HII recombination*/
			case 1: /* from a HeII recombination*/
			case 2: /* from a HeII recombination*/
			{
			const intptr_t ipar = r[i].ipar;
			for(d = 0; d < 3; d++) {
				r[i].s[d] = psys.pos[ipar][d];
			}
			const double T = psys_T(ipar);
			const double logT = log10(T);
			if(r[i].type == 0) {
				r[i].Nph = psys.yGrecHII[ipar] * psys_NH(ipar);
				r[i].freq = lte_gen_freq(LTE_FREQ_HI, logT);
			} else if(r[i].type == 1) {
				r[i].Nph = psys.yGrecHeII[ipar] * psys_NHe(ipar);
				r[i].freq = lte_gen_freq(LTE_FREQ_HEI, logT);
			} else if(r[i].type == 2) {
				r[i].Nph = psys.yGrecHeIII[ipar] * psys_NHe(ipar);
				r[i].freq = lte_gen_freq(LTE_FREQ_HEII, logT);
			}
			if(r[i].Nph < 0) ERROR("Nph =%g < 0, ipar= %ld", r[i].Nph, ipar);
			r[i].length = max_rec_length;
			gsl_ran_dir_3d(RNG, &dx, &dy, &dz);
			r[i].dir[0] = dx;
			r[i].dir[1] = dy;
			r[i].dir[2] = dz;
			break;
			}
			default:
			ERROR("never reach here");
		}
	}
	int identical_count = -10000;
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
			case -1: /* from a src */
				stat.src_ray_count.subtotal++;
				stat.src_photon_count.subtotal += r[i].Nph;
				psys.srcs[r[i].isrc].lastemit = psys.tick;
			break;
			case 0: /* from a recombination, aka particle */
			case 1: /* from a recombination, aka particle */
			case 2: /* from a recombination, aka particle */
			{
				const intptr_t ipar = r[i].ipar;
				stat.rec_ray_count.subtotal++;
				stat.rec_photon_count.subtotal += r[i].Nph;

				stat.rec_photon_count_sum -= r[i].Nph;
				if(r[i].Nph < 0.0) {
					WARNING("r[%ld].Nph = %g < 0.0", i, r[i].Nph);
				}
				if(r[i].type == 0) {
					stat.rec_photon_count.subtotalHII += r[i].Nph;
					psys.yGrecHII[ipar] -= r[i].Nph / psys_NH(ipar);
					if(psys.yGrecHII[ipar] > 1.0) {
						ERROR("ygrecHII[%ld] =%g > 0", ipar, psys.yGrecHII[ipar]);
					}
					if(psys.yGrecHII[ipar] < 0) {
						psys.yGrecHII[ipar] = 0;
					}
				} else if(r[i].type == 1) {
					stat.rec_photon_count.subtotalHeII += r[i].Nph;
					psys.yGrecHeII[ipar] -= r[i].Nph / psys_NHe(ipar);
					if(psys.yGrecHeII[ipar] < 0) {
						psys.yGrecHeII[ipar] = 0;
					}
				} else if(r[i].type == 2) {
					stat.rec_photon_count.subtotalHeIII += r[i].Nph;
					psys.yGrecHeIII[ipar] -= r[i].Nph / psys_NHe(ipar);
					if(psys.yGrecHeIII[ipar] < 0) {
						psys.yGrecHeIII[ipar] = 0;
					}
				}
			}
			break;
			default: ERROR("never each here");
		}
	}

	if(src_ran != NULL) gsl_ran_discrete_free(src_ran);
	stat.emit_time += omp_get_wtime() - t0;
}

static void trace() {
	intptr_t i;
	const double t0 = omp_get_wtime();
	#pragma omp parallel for private(i)
	for(i = 0; i < r_length; i++) {
		r[i].x_length = rt_trace(r[i].s, r[i].dir, r[i].length, 
			&r[i].x, &r[i].x_size);

		qsort(r[i].x, r[i].x_length, sizeof(Xtype), x_t_compare);
	}
	stat.raytrace_time += omp_get_wtime() - t0;
}

static void deposit(){
	const double t0 = omp_get_wtime();
	intptr_t i;

	const double scaling_fac = CFG_COMOVING?1/(psys.epoch->redshift + 1):1.0;
	const double scaling_fac2_inv = 1.0 / (scaling_fac * scaling_fac);
	const double scaling_fac3_inv = 1.0 / (scaling_fac * scaling_fac * scaling_fac);
	double spinlock_time = 0.0;
	size_t saturated_deposit_count = 0;
	double first_ionization_HI = 0;
	double first_ionization_HeI = 0;
	double first_ionization_HeII = 0;
	double secondary_ionization_HI = 0;
	double secondary_ionization_HeI = 0;
	size_t total_deposit_count = 0;
	size_t disordered_count = 0;
	#pragma omp parallel private(i) reduction(+: secondary_ionization_HI, secondary_ionization_HeI, first_ionization_HI, first_ionization_HeI, first_ionization_HeII, spinlock_time, total_deposit_count, saturated_deposit_count, disordered_count) 
	#pragma omp for schedule(dynamic, 1)
	for(i = 0; i < r_length; i++) {
		if(r[i].x_length != 0 && r[i].Nph != 0)
		{
		double Tau = 0.0;
		double TM = r[i].Nph; /*transmission*/
		if(CFG_NO_PHOTON && r[i].type < 0) TM = 0.0;
		intptr_t j;
		const double sigmaHI = xs_get(XS_HI, r[i].freq);
		const double sigmaHeI = xs_get(XS_HEI, r[i].freq);
		const double sigmaHeII = xs_get(XS_HEII, r[i].freq);
		const double dEHI = fdim(r[i].freq, C_HI_ENERGY);
		const double dEHeI = fdim(r[i].freq, C_HEI_ENERGY);
		const double dEHeII = fdim(r[i].freq, C_HEII_ENERGY);
		const double heat_factor_HI = C_H_PER_MASS * dEHI;
		const double heat_factor_HeI = C_HE_PER_MASS * dEHeI;
		const double heat_factor_HeII = C_HE_PER_MASS * dEHeII;

		total_deposit_count += r[i].x_length;

		for(j = 0; j < r[i].x_length; j++) {

			const intptr_t ipar = r[i].x[j].ipar;
			const double b = r[i].x[j].b;
			const double sml = psys.sml[ipar];
			const double sml_inv = 1.0 / sml;
			double NH = psys_NH(ipar);
			double NHe = psys_NHe(ipar);
			double rho = psys.rho[ipar] * scaling_fac3_inv;
			if (CFG_ENABLE_EOS) {
				const double x = eos_get_cloud_fraction(rho);
				rho *= (1 - x);
				NH *= (1 - x);
				NHe *= (1 - x);
			}
			const double NH_inv = 1.0 / NH;
			const double NHe_inv = 1.0 / NHe;
			const double xHI = psys_xHI(ipar);
			const double xHeI = psys_xHeI(ipar);
			const double xHeII = psys_xHeII(ipar);
			const double kernel = sph_depth(b * sml_inv) * (sml_inv * sml_inv) * scaling_fac2_inv;


			double deltaHI = 0.0;
			double deltaHeI = 0.0;
			double deltaHeII = 0.0;
			double tausum = 0.0;
			double tauHI = 0.0;
			double tauHeI = 0.0;
			double tauHeII = 0.0;

			const double NHI = xHI * NH;

			const double NHIcd = kernel * NHI;

			tauHI = sigmaHI * NHIcd;
			tausum += tauHI;

			if(!CFG_H_ONLY && NHe != 0.0) {
				const double NHeI = xHeI * NHe;
				const double NHeIcd = kernel * NHeI;
				tauHeI = sigmaHeI * NHeIcd;
				tausum += tauHeI;

				const double NHeII = xHeII * NHe;
				const double NHeIIcd = kernel * NHeII;
				tauHeII = sigmaHeII * NHeIIcd;
				tausum += tauHeII;
			}

			if(tausum == 0.0) {
				continue;
			}

			double absorb_est = -TM * expm1(-tausum);
			double absorb = 0;
			double absorbHI = tauHI / tausum * absorb_est;
			if(absorbHI > NHI) {
				absorbHI = NHI;
				saturated_deposit_count ++;
			}
			deltaHI = absorbHI * NH_inv;
			absorb = absorbHI;
			if(!CFG_H_ONLY && NHe != 0.0) {
				const double NHeI = xHeI * NHe;
				const double NHeII = xHeII * NHe;
				double absorbHeI = tauHeI / tausum * absorb_est;
				if(absorbHeI > NHeI) {
					absorbHeI = NHeI;
					saturated_deposit_count ++;
				}
				deltaHeI = absorbHeI * NHe_inv;
				double absorbHeII = tauHeII / tausum * absorb_est;
				if(absorbHeII > NHeII) {
					absorbHeII = NHeII;
					saturated_deposit_count ++;
				}
				deltaHeII = absorbHeII * NHe_inv;
				absorb += (absorbHeI + absorbHeII);
			}

			if(deltaHI < 0) {
				WARNING("absorbHI = %g NH_inv=%g NHI=%g\n tauHI=%g absorb_est=%g", absorbHI, NH_inv, NHI, tauHI, absorb_est);
				WARNING("tausum = %g, TM=%g nph %g", tausum, TM, r[i].Nph);
				ERROR("deltaHI < 0");
			}
			if(deltaHeI < 0) {
				ERROR("deltaHeI < 0");
			}
			if(deltaHeII < 0) {
				ERROR("deltaHeII < 0");
			}
			if(psys.heat[ipar] < 0.0) {
				ERROR("heat < 0");
			}

			double x = psys_ye(ipar) * NH / (NH + NHe);
			if(x > 1.0) x = 1.0;
			/* if x <= 0.0 secondary ionization has no effect */
			if(CFG_SECONDARY_IONIZATION && x > 0.0) {
#pragma omp atomic
				psys.yGdepHI[ipar] += deltaHI;
				first_ionization_HI += deltaHI * NH;
				double log10x = log10(x);
				double secion = deltaHI * dEHI / C_HI_ENERGY * secion_get(SECION_PHI_HI, dEHI, log10x)
					+ deltaHeI * dEHeI / C_HI_ENERGY * secion_get(SECION_PHI_HI, dEHeI, log10x);
					+ deltaHeII * dEHeII / C_HI_ENERGY * secion_get(SECION_PHI_HI, dEHeII, log10x);
#pragma omp atomic
				psys.yGdepHI[ipar] += secion;
				secondary_ionization_HI += secion * NH;
#pragma omp atomic
				psys.heat[ipar] += heat_factor_HI * deltaHI * secion_get(SECION_EH, dEHI, log10x);
				if(!CFG_H_ONLY && NHe != 0.0) {
					double secion = 
						deltaHI * dEHI / C_HEI_ENERGY * secion_get(SECION_PHI_HEI, dEHI, log10x)
						+ deltaHeI * dEHeI / C_HEI_ENERGY * secion_get(SECION_PHI_HEI, dEHeI, log10x);
						+ deltaHeII * dEHeII / C_HEI_ENERGY * secion_get(SECION_PHI_HEI, dEHeII, log10x);
#pragma omp atomic
					psys.yGdepHeI[ipar] += deltaHeI + secion;
					secondary_ionization_HeI += secion * NHe;
					first_ionization_HeI += deltaHeI * NHe;
#pragma omp atomic
					psys.heat[ipar] += heat_factor_HeI * deltaHeI * secion_get(SECION_EH, dEHeI, log10x);
						+ heat_factor_HeII * deltaHeII * secion_get(SECION_EH, dEHeII, log10x);
#pragma omp atomic
					psys.yGdepHeII[ipar] += deltaHeII;
					first_ionization_HeII += deltaHeII * NHe;
				}
			} else {
#pragma omp atomic
				psys.yGdepHI[ipar] += deltaHI;
#pragma omp atomic
				psys.heat[ipar] += heat_factor_HI * deltaHI;
				if(!CFG_H_ONLY && NHe != 0.0) {
#pragma omp atomic
					psys.yGdepHeI[ipar] += deltaHeI;
#pragma omp atomic
					psys.yGdepHeII[ipar] += deltaHeII;
#pragma omp atomic
					psys.heat[ipar] += deltaHeI * heat_factor_HeI;
#pragma omp atomic
					psys.heat[ipar] += deltaHeII * heat_factor_HeII;
				}
			}
			if(psys.heat[ipar] < 0.0) {
				ERROR("heat < 0");
			}
#pragma omp atomic
			psys.hits[ipar]++;


			maybe_write_particle(ipar, stat.hitlogfile, 
						"%lu %ld %g %g %g %g %g %g %g %g %g\n",
						psys.tick, ipar, xHI, b, sml, kernel, sigmaHI, sigmaHeI, sigmaHeII, NHI, absorb);

			TM -= absorb;
			if(TM / r[i].Nph < -1e-1) {
				ERROR("TM < 0 %g %g", TM, r[i].Nph);
			}
			/* cut off at around 10^-10 */
			if(TM < 0.0 || TM / r[i].Nph < 1e-10) {
				/* point to the next intersection so that a few lines later we can use j - 1*/
				j++;
				break;
			}
		}

//		ARRAY_RESIZE(r[i].x, Xtype, j);
//		r[i].length = r[i].x[j - 1].d;
		stat.lost_photon_count_sum += TM;
		}
	}
	stat.spinlock_time += spinlock_time;
	stat.deposit_time += (omp_get_wtime() - t0);
	stat.total_deposit_count += total_deposit_count;
	stat.saturated_deposit_count += saturated_deposit_count;
	stat.secondary_ionization.HI += secondary_ionization_HI;
	stat.secondary_ionization.HeI += secondary_ionization_HeI;
	stat.first_ionization.HI += first_ionization_HI;
	stat.first_ionization.HeI += first_ionization_HeI;
	stat.first_ionization.HeII += first_ionization_HeII;
	stat.disordered_count += disordered_count;
}

static void merge_pars() {
	const double t0 = omp_get_wtime();
	intptr_t i;
	/* now merge the list */
	size_t x_length_est = 0;
	for(i = 0; i < r_length; i++) {
		x_length_est += r[i].x_length;
	}

	ARRAY_ENSURE(active_set, intptr_t, x_length_est);
	ARRAY_CLEAR(active_set);

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
			if(bitmask_test(active, ipar)) {
				* (ARRAY_APPEND(active_set, intptr_t)) = r[i].x[j].ipar;
				bitmask_clear_unsafe(active, ipar);
			}
		}
		if(i <= psys.epoch->nray) {
			x_src_length = active_set_length;
		}
	}
//	qsort(active_set, x_src_length, sizeof(intptr_t), intptr_compare);
//	qsort(active_set + x_src_length, active_set_length - x_src_length, sizeof(intptr_t), intptr_compare);
	stat.merge_time += omp_get_wtime() - t0;
}

static void update_pars() {
	size_t d1 = 0, d2 = 0;
	double increase_recomb = 0;
	const double scaling_fac = CFG_COMOVING?1/(psys.epoch->redshift + 1.0):1.0;
	const double scaling_fac3_inv = 1.0/(scaling_fac * scaling_fac * scaling_fac);
	const double t0 = omp_get_wtime();
	#pragma omp parallel reduction(+: d1, d2, increase_recomb)
	{
	void * control = step_evolve_prepare();
	int NTH = omp_get_num_threads();
	int iTH = omp_get_thread_num();
	intptr_t j;
	for(j = 0; j < active_set_length + NTH; j+=NTH) {
		if(j + iTH >= active_set_length) continue;
		const intptr_t ipar = active_set[j + iTH];
		const double xHI = psys_xHI(ipar);
		const double xHII = psys_xHII(ipar);
		double rho = psys.rho[ipar] * scaling_fac3_inv;
		double NH = psys_NH(ipar);
		double NHe = psys_NHe(ipar);
		
		Step step = {0};

		step.ipar = ipar;
		step.time = (psys.tick - psys.lasthit[ipar]) * psys.tick_time;

		if (CFG_ENABLE_EOS) {
		/* EOS really only make sense with adiabatic simulations, because we can't model
         * the heating to the hot ambient in a meaningful way: it's internal energy only depends
         * on the density in the two-phase model */
			const double x = eos_get_cloud_fraction(rho);
			NH *= (1 - x);
			NHe *= (1 - x);
			if(CFG_FAKE_TEMPERATURE>0) {
				step.T = CFG_FAKE_TEMPERATURE;
			} else {
				if( x > 0.0) {
					step.T = ieye2T(eos_get_egyhot(rho), psys_ye(ipar));
				} else {
					step.T = psys_T(ipar);
				}
			}
		} else {
			if(CFG_FAKE_TEMPERATURE>0) {
				step.T = CFG_FAKE_TEMPERATURE;
			} else {
				step.T = psys_T(ipar);
			}
		}
		step.yGdepHI = psys.yGdepHI[ipar];
		step.yGdepHeI = psys.yGdepHeI[ipar];
		step.yGdepHeII = psys.yGdepHeII[ipar];
		step.heat = psys.heat[ipar];

		step.yGrecHII = psys.yGrecHII[ipar];
		step.yGrecHeII = psys.yGrecHeII[ipar];
		step.yGrecHeIII = psys.yGrecHeIII[ipar];

		step.lambdaH = psys.lambdaH[ipar];
		step.lambdaHeI = psys.lambdaHeI[ipar];
		step.lambdaHeII = psys.lambdaHeII[ipar];

		step.yeMET = psys.yeMET[ipar];
		step.nH = C_H_PER_MASS * rho;
		step.ie = psys.ie[ipar];

		maybe_write_particle(ipar, stat.parlogfile, "%lu %lu %g %g %g %g %g %g %g\n",
				psys.tick, ipar, step.time, step.T, step.lambdaH, step.nH, step.yGdepHI, step.ie, step.heat);

		if(!step_evolve(control, &step)) {
			const double xHI = lambdaH_to_xHI(step.lambdaH);
			const double xHII = lambdaH_to_xHII(step.lambdaH);
			WARNING("evolve failed: "
				"time = %g\n"
				"step.time = %g; \n"
				"step.T = %g; \n"
				"step.lambdaH=%g; \n"
				"step.yeMET=%g;\n"
				"step.nH=%g;\n"
				"step.yGdep=%g;\n"
				"step.ie=%g;\n"
				"step.heat=%g;\n"
				"xHI=%g\n"
				"xHII=%g\n",
				step.time/U_MYR, step.time, step.T, 
				step.lambdaH, 
				step.yeMET, step.nH, step.yGdepHI, step.ie, step.heat,
				xHI, xHII
				);
			d1++;
		} else {
			if(step.refined){
				stat.fast_recombination_count++;
			}
			increase_recomb += (step.yGrecHII - psys.yGrecHII[ipar])* NH;
			increase_recomb += (step.yGrecHeII - psys.yGrecHeII[ipar])* NHe;
			increase_recomb += (step.yGrecHeIII - psys.yGrecHeIII[ipar])* NHe;

			psys.yGdepHI[ipar] *= step.step_remain;
			psys.yGdepHeI[ipar] *= step.step_remain;
			psys.yGdepHeII[ipar] *= step.step_remain;

			psys.yGrecHII[ipar] = step.yGrecHII;
			psys.yGrecHeII[ipar] = step.yGrecHeII;
			psys.yGrecHeIII[ipar] = step.yGrecHeIII;

			psys.lambdaH[ipar] = fmax(0, fmin(1, step.lambdaH));
			psys.lambdaHeI[ipar] = fmax(0, fmin(1, step.lambdaHeI));
			psys.lambdaHeII[ipar] = fmax(0, fmin(1 - psys.lambdaHeI[ipar], step.lambdaHeII));

			if(!CFG_ADIABATIC)
				psys.ie[ipar] = step.ie;

			psys.lasthit[ipar] = psys.tick;

			if(CFG_ISOTHERMAL) {
				psys.ie[ipar] = Tye2ie(step.T, psys_ye(ipar));
			}

			psys.heat[ipar] = 0.0;
		}
		d2++;
	}
	step_evolve_free(control);
	}
	stat.gsl_error_count += d1;
	stat.evolve_count += d2;
	stat.rec_photon_count_sum += increase_recomb;
	stat.update_time += omp_get_wtime() - t0;
}

static int intptr_compare(const void * p1, const void * p2) {
	const intptr_t d1 = *((const intptr_t * )p1);
	const intptr_t d2 = *((const intptr_t * )p2);
    return (d1 > d2) - (d2 > d1);
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
