#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <messages.h>
#include <array.h>
#include "config.h"
#include "reader.h"
#include "psystem.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

static void solve_u_v(double d[3], double u[3], double v[3]);
static void psystem_read_source_file(const char * filename);
static void psystem_read_source_file2(const char * filename);

void psystem_read_source() {
	if(psys.epoch->sources) {
		const double scaling_fac = CFG_COMOVING?1/(psys.epoch->redshift + 1):1.0;
		if(psys.srcs) free(psys.srcs);
		psys.srcs = NULL;
		psys.nsrcs = 0;
		int i = 0;
		for(i = 0; i < psys.epoch->nsources; i++) {
			psystem_read_source_file2(psys.epoch->sources[i]);
		}
		intptr_t isrc;
		for(isrc = 0; isrc < psys.nsrcs; isrc++) {
			psys.srcs[isrc].ray_length_hint = C_SPEED_LIGHT / scaling_fac * psys.tick_time;
		}
		double Lmin;
		double Lmax;
		intptr_t imin = -1, imax = -1;
		for(isrc = 0; isrc < psys.nsrcs; isrc++) {
			if(isrc == 0 || psys_Ngamma_dot(isrc) > Lmax) {
				Lmax = psys_Ngamma_dot(isrc);
				imax = isrc;
			}
			if(isrc == 0 || psys_Ngamma_dot(isrc) < Lmin) {
				Lmin = psys_Ngamma_dot(isrc);
				imin = isrc;
			}
		}
		MESSAGE("SOURCES: %lu, max=%g at (%g %g %g), min=%g at (%g %g %g)",
			psys.nsrcs, Lmax * U_SEC, 
			psys.srcs[imax].pos[0],
			psys.srcs[imax].pos[1],
			psys.srcs[imax].pos[2],
			Lmin * U_SEC,
			psys.srcs[imin].pos[0],
			psys.srcs[imin].pos[1],
			psys.srcs[imin].pos[2]);
	}
}

static void advance_cursor(intptr_t i) {
	while( psys.srcs[i].cursor < psys.srcs[i].ticks_length - 1
	    && psys.srcs[i].ticks[psys.srcs[i].cursor + 1] <= psys.tick
	) {
		psys.srcs[i].cursor++;
		MESSAGE("advance source %ld cursor to %ld, tick = %ld",
			i, psys.srcs[i].cursor, psys.tick);
	}
}
void psystem_weight_srcs(double weights[]) {
	intptr_t i;
	for(i = 0; i < psys.nsrcs; i++) {
		advance_cursor(i);
		/* treat two types the same essentially because they are both Ngamma_sec*/
		if(psys.srcs[i].type == PSYS_SRC_POINT) {
			weights[i] = psys_Ngamma_dot(i) * (psys.tick - psys.srcs[i].lastemit) * psys.tick_time;
		} else if(psys.srcs[i].type == PSYS_SRC_PLANE) {
			weights[i] = psys_Ngamma_dot(i) * (psys.tick - psys.srcs[i].lastemit) * psys.tick_time;
		}
	}
}

static void solve_u_v(double d[3], double u[3], double v[3]) {
	double data[9] = {
		d[0] * d[0], d[1] * d[0], d[2] * d[0],
		d[0] * d[1], d[1] * d[1], d[2] * d[1],
		d[0] * d[2], d[1] * d[2], d[2] * d[2],
	};
	gsl_matrix_view m = gsl_matrix_view_array(data, 3, 3);
	gsl_vector * eval = gsl_vector_alloc(3);
	gsl_matrix * evac = gsl_matrix_alloc(3, 3);
	gsl_eigen_symmv_workspace * work = gsl_eigen_symmv_alloc(3);
	gsl_eigen_symmv(&m.matrix, eval, evac, work);
	gsl_eigen_symmv_sort(eval, evac, GSL_EIGEN_SORT_ABS_ASC);
	gsl_vector_view u_view = gsl_matrix_column(evac, 0);
	gsl_vector_view v_view = gsl_matrix_column(evac, 1);
	int i;
	for(i = 0; i < 3; i++) {
		u[i] = gsl_vector_get(&u_view.vector, i);
		v[i] = gsl_vector_get(&v_view.vector, i);
	}
	gsl_matrix_free(evac);
	gsl_vector_free(eval);
	gsl_eigen_symmv_free(work);
}

static void psystem_read_source_file(const char * filename) {
	FILE * f = fopen(filename, "r");
	if(f == NULL) {
		ERROR("failed to open %s", filename);
	}
	MESSAGE("reading a version 1 source file");
	int NR = 0;
	int NF = 0;
	char * line = NULL;
	size_t len = 0;
	intptr_t isrc = 0;
	int stage = 0;
	double x, y, z, dx, dy, dz, radius, L;
	char spec[128];
	char type[128];

	size_t nsrcs_file;
	while(0 <= getline(&line, &len, f)) {
		if(line[0] == '#' || line[0] == '\n') {
			NR++;
			continue;
		}
		switch(stage) {
		case 0:
			if(1 != sscanf(line, "%ld", &nsrcs_file)) {
				ERROR("%s format error at %d", filename, NR);
			} else {
				isrc = psys.nsrcs;
				psys.nsrcs += nsrcs_file;
				psys.srcs = realloc(psys.srcs, sizeof(Source) * (psys.nsrcs));
				stage ++;
			}
			break;
		case 1:
			NF = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %80s %80s %lf",
				&x, &y, &z, &dx, &dy, &dz, 
				&L, spec, type, &radius); 
			if(NF != 9 && NF != 10) {
				ERROR("%s format error at %d", filename, NR);
			}
			psys.srcs[isrc].pos[0] = x;
			psys.srcs[isrc].pos[1] = y;
			psys.srcs[isrc].pos[2] = z;
			psys.srcs[isrc].dir[0] = dx;
			psys.srcs[isrc].dir[1] = dy;
			psys.srcs[isrc].dir[2] = dz;
			psys.srcs[isrc].specid = spec_get(spec);
			psys.srcs[isrc].ticks = NULL;
			psys.srcs[isrc].Ngamma_dots = NULL;
			psys.srcs[isrc].lastemit = 0;
			ARRAY_RESIZE(psys.srcs[isrc].ticks, intptr_t, 1);
			ARRAY_RESIZE(psys.srcs[isrc].Ngamma_dots, double, 1);
			psys.srcs[isrc].cursor = 0;
			if(!strcmp(type, "plane")) {
				if(NF != 10) {
					ERROR("%s format error at %d, needs 10 fields", filename, NR);
				}
				psys.srcs[isrc].type = PSYS_SRC_PLANE;
				psys.srcs[isrc].radius = radius;
				psys.srcs[isrc].ticks[0] = 0;
				psys.srcs[isrc].Ngamma_dots[0] = L * M_PI * radius * radius / (U_CM * U_CM * U_SEC);
				solve_u_v(psys.srcs[isrc].dir, psys.srcs[isrc].a, psys.srcs[isrc].b);
			} else {
				if(L < 0) {
					L = spec_N_from_lum(psys.srcs[isrc].specid, -L * C_SOLAR_LUM) * U_SEC / 1e50;
				}
				psys.srcs[isrc].type = PSYS_SRC_POINT;
				psys.srcs[isrc].ticks[0] = 0;
				psys.srcs[isrc].Ngamma_dots[0] = L * 1e50 / U_SEC;
			}
			psys.srcs[isrc].Ngamma_dots[0] *= CFG_BOOST_SOURCE_FACTOR; /* boost the photon luminosity according to cfg file */
			isrc ++;
			if(isrc == psys.nsrcs) {
				stage ++;
			}
			break;
		case 2:
			break;
		}
		NR++;
		if(stage == 2) break;
	}
	if(isrc != psys.nsrcs) {
		ERROR("%s expecting %lu, but has %lu sources, %lu", filename, psys.nsrcs, isrc, NR);
	}
	free(line);
	fclose(f);
}
static void psystem_read_source_file2(const char * filename) {
	FILE * f = fopen(filename, "r");
	if(f == NULL) {
		ERROR("failed to open %s", filename);
	}
	int NR = 0;
	int NF = 0;
	char * line = NULL;
	size_t len = 0;
	intptr_t isrc = 0;
	int stage = 0;
	size_t Nsamples;
	intptr_t isample;
	double x, y, z, dx, dy, dz, radius, L, time;
	char spec[128];
	char type[128];

	char timefmtstr[128];
	char locationfmtstr[128] = "absolute";
	int timefmt = 0;
	int locationfmt = 0;
	size_t nsrcs_file;

	while(0 <= getline(&line, &len, f)) {
		if(line[0] == '#' || line[0] == '\n') {
			NR++;
			continue;
		}
		switch(stage) {
		case 0:
			if(1 == sscanf(line, "%ld %80s %80s", &nsrcs_file, timefmtstr, locationfmtstr)) {
				fclose(f);
				free(line);
				/* maybe this is a verion 1 file? */
				psystem_read_source_file(filename);
				return;
			} else {
				if(!strcmp(timefmtstr, "time_q")) 
					timefmt = 0;
				else if(!strcmp(timefmtstr, "a")) 
					timefmt = 0;
				else if(!strcmp(timefmtstr, "redshift")) 
					timefmt = 1;
				else if(!strcmp(timefmtstr, "z")) 
					timefmt = 1;
				else if(!strcmp(timefmtstr, "myr")) 
					timefmt = 2;
				else if(!strcmp(timefmtstr, "tick")) 
					timefmt = 3;
				else if(!strcmp(locationfmtstr, "absolute")) {
					locationfmt = 0;
				} else {
					ERROR("time fmt str %s unknown", timefmtstr);
				}
				if(!strcmp(locationfmtstr, "relative")) {
					locationfmt = 1;
				}
				isrc = psys.nsrcs;
				psys.nsrcs += nsrcs_file;
				psys.srcs = realloc(psys.srcs, sizeof(Source) * psys.nsrcs);
				stage ++;
			}
			break;
		case 1:
			NF = sscanf(line, "%ld %80s %80s %lf", &Nsamples, spec, type, &radius);
			if(NF == 3 && !strcmp(type, "plane")) {
				ERROR("%s format error at %d, needs radius for plane source", filename, NR);
			}
			psys.srcs[isrc].ticks = NULL;
			psys.srcs[isrc].Ngamma_dots = NULL;
			psys.srcs[isrc].lastemit = 0;
			ARRAY_RESIZE(psys.srcs[isrc].ticks, intptr_t, Nsamples);
			ARRAY_RESIZE(psys.srcs[isrc].Ngamma_dots, double, Nsamples);
			psys.srcs[isrc].cursor = 0;
			psys.srcs[isrc].specid = spec_get(spec);
			if(!strcmp(type, "plane")) {
				psys.srcs[isrc].type = PSYS_SRC_PLANE;
				psys.srcs[isrc].radius = radius;
			} else {
				psys.srcs[isrc].type = PSYS_SRC_POINT;

			}
			stage++;
			isample = 0;
			break;
		case 2:
			NF = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf",
				&time, &x, &y, &z, &dx, &dy, &dz, &L ); 
			if(NF != 8) {
				ERROR("%s format error at %d", filename, NR);
			}
			intptr_t tick;

			switch(timefmt) {
				case 0:
					time = z2t(1 / (1 + time));
					tick = (time * U_MYR - psys.epoch->age) / psys.tick_time;
				case 1:
					time = z2t(time);
					tick = (time * U_MYR - psys.epoch->age) / psys.tick_time;
				break;
				case 2:
					tick = (time * U_MYR) / psys.tick_time;
				break;
				case 3:
					tick = time;
				break;
			}
			psys.srcs[isrc].ticks[isample] = tick;

			switch(locationfmt) {
				case 0:
				break;
				case 1:
				x = psys.boxsize * x;
				y = psys.boxsize * y;
				z = psys.boxsize * z;
				break;
			}

			if(psys.srcs[isrc].type == PSYS_SRC_PLANE) {
				psys.srcs[isrc].Ngamma_dots[isample] = L * M_PI * radius * radius / (U_CM * U_CM * U_SEC);
				solve_u_v(psys.srcs[isrc].dir, psys.srcs[isrc].a, psys.srcs[isrc].b);
			}
			if(psys.srcs[isrc].type == PSYS_SRC_POINT) {
				if(L < 0) {
					L = spec_N_from_lum(psys.srcs[isrc].specid, -L * C_SOLAR_LUM) * U_SEC / 1e50;
				}
				psys.srcs[isrc].Ngamma_dots[isample] = L * 1e50 / U_SEC;
			}
			psys.srcs[isrc].Ngamma_dots[isample] *= CFG_BOOST_SOURCE_FACTOR; /* boost the photon luminosity according to cfg file */

			psys.srcs[isrc].pos[0] = x;
			psys.srcs[isrc].pos[1] = y;
			psys.srcs[isrc].pos[2] = z;
			psys.srcs[isrc].dir[0] = dx;
			psys.srcs[isrc].dir[1] = dy;
			psys.srcs[isrc].dir[2] = dz;
			isample ++;
			if(isample == Nsamples) {
				if(isrc == psys.nsrcs) {
					stage ++;
				} else {
					stage --;
					isrc++;
				}
			}
			break;
		case 3:
			break;
		}
		NR++;
		if(stage == 3) break;
	}
	if(isrc != psys.nsrcs) {
		ERROR("%s expecting %lu, but has %lu sources, %lu", filename, psys.nsrcs, isrc, NR);
	}
	free(line);
	fclose(f);

}
