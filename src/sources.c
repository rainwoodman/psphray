#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <messages.h>

#include "config.h"
#include "reader.h"
#include "psystem.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

static void solve_u_v(double d[3], double u[3], double v[3]);
static void psystem_read_source_file(char * filename);

void psystem_read_source() {
	if(psys.epoch->source) {
		const double scaling_fac = CFG_COMOVING?1/(psys.epoch->redshift + 1):1.0;
		if(psys.srcs) free(psys.srcs);
		psystem_read_source_file(psys.epoch->source);
		intptr_t isrc;
		for(isrc = 0; isrc < psys.nsrcs; isrc++) {
			psys.srcs[isrc].ray_length_hint = C_SPEED_LIGHT / scaling_fac * psys.tick_time;
		}
		double Lmin;
		double Lmax;
		intptr_t imin = -1, imax = -1;
		for(isrc = 0; isrc < psys.nsrcs; isrc++) {
			if(isrc == 0 || psys.srcs[isrc].Ngamma_dot > Lmax) {
				Lmax = psys.srcs[isrc].Ngamma_dot ;
				imax = isrc;
			}
			if(isrc == 0 || psys.srcs[isrc].Ngamma_dot < Lmin) {
				Lmin = psys.srcs[isrc].Ngamma_dot;
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

void psystem_get_source_weights(double weights[]) {
	intptr_t i;
	for(i = 0; i < psys.nsrcs; i++) {
		/* treat two types the same essentially because they are both Ngamma_sec*/
		if(psys.srcs[i].type == PSYS_SRC_POINT) {
			weights[i] = psys.srcs[i].Ngamma_dot * (psys.tick - psys.srcs[i].lastemit) * psys.tick_time;
		} else if(psys.srcs[i].type == PSYS_SRC_PLANE) {
			weights[i] = psys.srcs[i].Ngamma_dot * (psys.tick - psys.srcs[i].lastemit) * psys.tick_time;
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

static void psystem_read_source_file(char * filename) {
	FILE * f = fopen(filename, "r");
	if(f == NULL) {
		ERROR("failed to open %s", psys.epoch->source);
	}
	int NR = 0;
	int NF = 0;
	char * line = NULL;
	size_t len = 0;
	intptr_t isrc = 0;
	int stage = 0;
	double x, y, z, dx, dy, dz, radius, L;
	char spec[128];
	char type[128];

	while(0 <= getline(&line, &len, f)) {
		if(line[0] == '#') {
			NR++;
			continue;
		}
		switch(stage) {
		case 0:
			if(1 != sscanf(line, "%ld", &psys.nsrcs)) {
				ERROR("%s format error at %d", psys.epoch->source, NR);
			} else {
				psys.srcs = calloc(sizeof(Source), psys.nsrcs);
				isrc = 0;
				stage ++;
			}
			break;
		case 1:
			NF = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %80s %80s %lf",
				&x, &y, &z, &dx, &dy, &dz, 
				&L, spec, type, &radius); 
			if(NF != 9 && NF != 10) {
				ERROR("%s format error at %d", psys.epoch->source, NR);
			}
			psys.srcs[isrc].pos[0] = x;
			psys.srcs[isrc].pos[1] = y;
			psys.srcs[isrc].pos[2] = z;
			psys.srcs[isrc].dir[0] = dx;
			psys.srcs[isrc].dir[1] = dy;
			psys.srcs[isrc].dir[2] = dz;
			psys.srcs[isrc].specid = spec_get(spec);
			if(!strcmp(type, "plane")) {
				if(NF != 10) {
					ERROR("%s format error at %d, needs 10 fields", psys.epoch->source, NR);
				}
				psys.srcs[isrc].type = PSYS_SRC_PLANE;
				psys.srcs[isrc].radius = radius;
				psys.srcs[isrc].Ngamma_dot = L * M_PI * radius * radius / (U_CM * U_CM / U_SEC);
				solve_u_v(psys.srcs[isrc].dir, psys.srcs[isrc].a, psys.srcs[isrc].b);
			} else {
				if(L < 0) {
					L = spec_N_from_lum(psys.srcs[isrc].specid, -L * C_SOLAR_LUM) * U_SEC / 1e50;
				}
				psys.srcs[isrc].type = PSYS_SRC_POINT;
				psys.srcs[isrc].Ngamma_dot = L * 1e50 / U_SEC;
			}
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
		ERROR("%s expecting %lu, but has %lu sources, %lu", psys.epoch->source, psys.nsrcs, isrc, NR);
	}
	free(line);
	fclose(f);
}
