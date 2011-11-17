#include <stdio.h>
#include <stdarg.h>

#include <messages.h>
#include <stdint.h>

#include "stat.h"
#include "config.h"
#include "psystem.h"

#include <omp.h>

stat_t stat = {0};

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

void stat_restart() {
	if(CFG_DUMP_HOTSPOTS) {
		stat.parlogfile = fopen_printf("parlogfile-%03d", "w", psys.epoch - EPOCHS);
		stat.hitlogfile = fopen_printf("hitlogfile-%03d", "w", psys.epoch - EPOCHS);
	}
}

void stat_subtotal() {
	stat.src_ray_count.total += stat.src_ray_count.subtotal;
	stat.rec_ray_count.total += stat.rec_ray_count.subtotal;
	stat.src_photon_count.total += stat.src_photon_count.subtotal;
	stat.rec_photon_count.total += stat.rec_photon_count.subtotal;

	MESSAGE("SRC(subtot): Ray %lu Photon %g ", 
		stat.src_photon_count.subtotal, stat.src_ray_count.subtotal);
	MESSAGE("REC(subtot): Ray %lu Photon %g ", 
		stat.rec_photon_count.subtotal, stat.rec_ray_count.subtotal);
	MESSAGE("REC(subtot): HII %g, HeII %g HeIII %g", 
		stat.rec_photon_count.subtotalHII, stat.rec_photon_count.subtotalHeII, stat.rec_photon_count.subtotalHeIII);
	MESSAGE("PH EMISSION: Src %g Rec %g",
		stat.src_photon_count.total, stat.rec_photon_count.total);
	MESSAGE("PH STORAGE : Lost %g Rec %g ", 
		stat.lost_photon_count_sum, stat.rec_photon_count_sum);
	MESSAGE("Deposit: saturated = %lu, disordered= %lu,  total=%lu", stat.saturated_deposit_count, stat.disordered_count, stat.total_deposit_count);
	MESSAGE("First:   %g %g %g", stat.first_ionization.HI, stat.first_ionization.HeI, stat.first_ionization.HeII);
	MESSAGE("Secondary:   %g %g", stat.secondary_ionization.HI, stat.secondary_ionization.HeI);
	MESSAGE("Evolve     : Error %lu FastRec %lu Total %lu", 
		stat.gsl_error_count, stat.fast_recombination_count, stat.evolve_count);
	MESSAGE("Time : SpinLock %g", stat.spinlock_time);

	stat.src_ray_count.subtotal = 0;
	stat.rec_ray_count.subtotal = 0;
	stat.src_photon_count.subtotal = 0;
	stat.rec_photon_count.subtotal = 0;
	stat.rec_photon_count.subtotalHII = 0;
	stat.rec_photon_count.subtotalHeII = 0;
	stat.rec_photon_count.subtotalHeIII = 0;
	MESSAGE("Real Time: emit %g raytrace %g deposit %g update %g total %g",
		stat.emit_time, stat.raytrace_time, stat.deposit_time, stat.update_time, stat.total_time);
	stat.tick_subtotal = 0;
}

void stat_stop() {
	if(CFG_DUMP_HOTSPOTS) {
		fclose(stat.parlogfile);
		fclose(stat.hitlogfile);
	}
}
