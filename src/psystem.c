#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <messages.h>

#include "config.h"
#include "reader.h"
#include "psystem.h"

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_heapsort.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#define IDHASHBITS 24
#define IDHASHMASK ((((size_t)1) << IDHASHBITS) - 1)


#define read_as(r, blk, ptype, dst, fromtype, totype) \
{ if(!strcmp(# fromtype, # totype)) { \
	reader_read(r, blk, ptype, dst); \
} else { \
	fromtype * buf = reader_alloc(r, blk, ptype); \
	reader_read(r, blk, ptype, buf); \
	size_t npar_file = reader_npar(r, ptype); \
	intptr_t i; \
	for(i = 0; i < npar_file; i++) { \
		(dst)[i] = buf[i]; \
	} \
	free(buf); \
} }

#define write_as(r, blk, ptype, dst, fromtype, totype) \
{ if(!strcmp(# fromtype , # totype)) { \
	reader_write(r, blk, ptype, dst); \
} else { \
	totype * buf = reader_alloc(r, blk, ptype); \
	size_t npar_file = reader_npar(r, ptype); \
	intptr_t i; \
	for(i = 0; i < npar_file; i++) { \
		buf[i] = (dst)[i]; \
	} \
	reader_write(r, blk, ptype, buf); \
	free(buf); \
} }

PSystem psys = {0};

static void idhash_build(uint64_t * id, size_t n) {
	psys.idhash.head = malloc(sizeof(intptr_t) * (((size_t)1) << IDHASHBITS));
	psys.idhash.next = malloc(sizeof(intptr_t) * n);
	intptr_t i;
	memset(psys.idhash.head, -1, sizeof(intptr_t) * (((size_t)1) << IDHASHBITS));
	memset(psys.idhash.next, -1, sizeof(intptr_t) * n);
	for(i = 0; i < n; i++) {
		intptr_t hash = id[i] & IDHASHMASK;
		psys.idhash.next[i] = psys.idhash.head[hash];
		psys.idhash.head[hash] = i;
	}
}
static int intptr_t_compare(const intptr_t * p1, const intptr_t * p2) {
	if(*p1 > *p2) return 1;
	if(*p1 < *p2) return -1;
	if(*p1 == *p2) return 0;
}

static void * permute (const size_t * p, void * data, const size_t ele_bytes, const size_t stride, const size_t n)
{
	void * out = malloc(stride * n);
	size_t i;
	#pragma omp parallel for private(i)
	for(i = 0; i < n; i++) {
		if(p[i] > n ) ERROR("p[i] out of range check peano key");
		memcpy((char*) out+ i * stride, (char*) data + p[i] * stride, ele_bytes);
	}
	free(data);
	return out;
}

static void hilbert_reorder() {
	extern intptr_t peano_hilbert_key(int x, int y, int z, int bits);
	size_t * perm = malloc(sizeof(size_t) * psys.npar);
	intptr_t * peanokeys = malloc(sizeof(intptr_t) * psys.npar);
	intptr_t ipar;
	float fac = 1 / (psys.boxsize * (1L << 20));
	#pragma omp parallel for private(ipar)
	for(ipar = 0; ipar < psys.npar; ipar++) {
		float * pos = psys.pos[ipar];
		peanokeys[ipar] = peano_hilbert_key(
			pos[0] * fac,
			pos[1] * fac,
			pos[2] * fac, 20);
	}
	gsl_heapsort_index(perm, peanokeys, psys.npar, sizeof(intptr_t), (void*)intptr_t_compare);

	psys.pos = permute(perm, psys.pos, 3 * sizeof(float), 3 * sizeof(float), psys.npar);
	psys.lambdaH = permute(perm, psys.lambdaH, sizeof(double), sizeof(double), psys.npar);
	psys.lambdaHeI = permute(perm, psys.lambdaHeI, sizeof(double), sizeof(double), psys.npar);
	psys.lambdaHeII = permute(perm, psys.lambdaHeII, sizeof(double), sizeof(double), psys.npar);
	psys.yeMET = permute(perm, psys.yeMET, sizeof(float), sizeof(float), psys.npar);
	psys.mass = permute(perm, psys.mass, sizeof(float), sizeof(float), psys.npar);
	psys.sml = permute(perm, psys.sml, sizeof(float), sizeof(float), psys.npar);
	psys.rho = permute(perm, psys.rho, sizeof(float), sizeof(float), psys.npar);
	psys.ie = permute(perm, psys.ie, sizeof(float), sizeof(float), psys.npar);
	psys.flag = permute(perm, psys.flag, sizeof(int8_t), sizeof(int8_t), psys.npar);
	int i;
	for(i = 0; i < 3; i++) {
		psys.yGrec[i] = permute(perm, psys.yGrec[i], sizeof(float), sizeof(float), psys.npar);
		psys.yGdep[i] = permute(perm, psys.yGdep[i], sizeof(float), sizeof(float), psys.npar);
	}
	psys.heat = permute(perm, psys.heat, sizeof(float), sizeof(float), psys.npar);
	psys.lasthit = permute(perm, psys.lasthit, sizeof(intptr_t), sizeof(intptr_t), psys.npar);
	psys.hits = permute(perm, psys.hits, sizeof(int), sizeof(int), psys.npar);
	psys.id = permute(perm, psys.id, sizeof(uint64_t), sizeof(uint64_t), psys.npar);

	free(peanokeys);
	free(perm);
}

static float dist(const float p1[3], const float p2[3]);
static void psystem_read_source();
static void psystem_read_epoch(ReaderConstants * c);
static void psystem_match_epoch(ReaderConstants * c);

static void setup_hotspots() {
	intptr_t ipar;
	#pragma omp parallel for private(ipar)
	for(ipar = 0; ipar < psys.npar; ipar++) {
		psys.flag[ipar] &= ~PF_HOTSPOT;
	}
	config_setting_t * hotspots = config_lookup(CFG, "psphray.hotspots");
	if(hotspots) {
		size_t n = config_setting_length(hotspots);
		intptr_t i;
		for(i = 0; i < n; i++) {
			config_setting_t * hotspot = config_setting_get_elem(hotspots, i);
			config_setting_t * center = config_setting_get_elem(hotspot, 0);
			config_setting_t * radius = config_setting_get_elem(hotspot, 1);
			int d;
			double c[3], r;
			for(d = 0; d < 3; d++) {
				config_setting_parse_units_elem(center, d, &c[d]);
			}
			r = config_setting_parse_units(radius);
			intptr_t ipar;
			#pragma omp parallel for private(ipar)
			for(ipar = 0; ipar < psys.npar; ipar++) {
				double dist = 0;
				int d;
				for(d = 0; d < 3; d++) {
					const double dx = (psys.pos[ipar][d] - c[d]);
					dist += dx * dx;
				}
				if(dist < r * r) {
					psys.flag[ipar] |= PF_HOTSPOT;
				}
			}
		}
	}
}
void psystem_switch_epoch(int i) {

	psys.epoch = &EPOCHS[i];

	psys.tick = 0;
	psys.tick_time = psys.epoch->duration / psys.epoch->nticks;

	Reader * r0 = reader_new(psys.epoch->format);
	char * fname0 = reader_make_filename(psys.epoch->snapshot, 0);
	reader_open(r0, fname0);
	free(fname0);
	ReaderConstants * c = reader_constants(r0);

	/* read stuff */

	config_lookup_float(CFG, "box.boxsize", &psys.boxsize);
	if(psys.boxsize == -1) {
		psys.boxsize = c->boxsize;
	}
	const char * boundary = NULL;
	config_lookup_string(CFG, "box.boundary", &boundary);
	if(!strcmp(boundary, "vaccum"))
		psys.periodic = 0;
	else if(!strcmp(boundary, "periodic"))
		psys.periodic = 1;
	else ERROR("boundary type %s unknown, be vaccum or periodic", boundary);
	MESSAGE("epoch %d, redshift %g, snapshot %s, format %s, source %s, ticks %u, ngas %lu, age %f [myr], duration %f [myr], boxsize = %g",
	i, 
	psys.epoch->redshift,
	psys.epoch->snapshot,
	psys.epoch->format,
	psys.epoch->source,
	psys.epoch->nticks,
	psys.epoch->ngas,
	units_format(psys.epoch->age, "myr"),
	units_format(psys.epoch->duration, "myr"),
	psys.boxsize);

	if(i == 0) {
		/* initial read */
		psystem_read_epoch(c);
	} else {
		/* matching with existing psystem */
		psystem_match_epoch(c);
	}

	setup_hotspots();

	hilbert_reorder();

	/*FIXME: calculate active particles*/
	MESSAGE("EPOCH active gas particles %ld/%ld", 0, c->Ntot[0]);

	intptr_t ipar;
	double mass = 0;
	for(ipar = 0; ipar < psys.npar; ipar++) {
		mass += psys.mass[ipar];
	}
	MESSAGE("EPOCH # of protons %g", C_HMF * mass / U_MPROTON);

	/* free r0 here, because we need c till now*/
	reader_destroy(r0);

	if(psys.epoch->source) {
		const double scaling_fac = CFG_COMOVING?1/(psys.epoch->redshift + 1):1.0;
		if(psys.srcs) free(psys.srcs);
		psystem_read_source();
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

static void psystem_read_epoch(ReaderConstants * c) {
	size_t ngas = psys.epoch->ngas;
	psys.npar = ngas;
	int fid;
	size_t nread = 0;
	psys.pos = calloc(sizeof(float) * 3, ngas);
	psys.mass = calloc(sizeof(float), ngas);
	psys.id = calloc(sizeof(uint64_t), ngas);
	psys.ie = calloc(sizeof(float), ngas);
	psys.sml = calloc(sizeof(float), ngas);
	psys.rho = calloc(sizeof(float), ngas);
	psys.lambdaH = calloc(sizeof(double), ngas);
	psys.lambdaHeI = calloc(sizeof(double), ngas);
	psys.lambdaHeII = calloc(sizeof(double), ngas);
	psys.yeMET = calloc(sizeof(float), ngas);
	int i;
	for(i = 0; i < 3; i++) {
		psys.yGrec[i] = calloc(sizeof(float), ngas);
		psys.yGdep[i] = calloc(sizeof(float), ngas);
	}
	psys.heat = calloc(sizeof(float), ngas);
	psys.lasthit = calloc(sizeof(intptr_t), ngas);
	psys.hits = calloc(sizeof(int), ngas);

	psys.flag = calloc(sizeof(int8_t), ngas);

	for(fid = 0; fid < c->Nfiles; fid ++) {
		Reader * r = reader_new(psys.epoch->format);
		char * fname = reader_make_filename(psys.epoch->snapshot, fid);
		reader_open(r, fname);
		free(fname);
		reader_read(r, "pos", 0, psys.pos[nread]);
		reader_read(r, "mass", 0, &psys.mass[nread]);
		reader_read(r, "sml", 0, &psys.sml[nread]);
		reader_read(r, "rho", 0, &psys.rho[nread]);
		reader_read(r, "ie", 0, &psys.ie[nread]);

		reader_read(r, "ye", 0, &psys.yeMET[nread]);
		float * xHI = reader_alloc(r, "xHI", 0);
		reader_read(r, "xHI", 0, xHI);

		intptr_t ipar;
		size_t npar_file = reader_npar(r, 0);

		for(ipar = 0; ipar < npar_file; ipar++) {
			const double xHII = 1.0 - xHI[ipar];
			psys_set_lambdaH(nread + ipar, xHI[ipar], xHII);

			const double yeMET = psys.yeMET[nread + ipar] - xHII;
			psys.yeMET[nread + ipar] = (yeMET < 0.0)?0.0:yeMET;
		}

		free(xHI);
		/* initialize He to fully neutral, for now */
		for(ipar = 0; ipar < npar_file; ipar++) {
			psys_set_lambdaHe(nread + ipar, 1.0, 0.0, 0.0);
		}

		if(reader_itemsize(r, "id") == 4) {
			unsigned int * id = reader_alloc(r, "id", 0);
			reader_read(r, "id", 0, id);
			intptr_t ipar;
			for(ipar = 0; ipar < npar_file; ipar++) {
				psys.id[nread + ipar] = id[ipar];
			}
			free(id);
		} else {
			reader_read(r, "id", 0, &psys.id[nread]);
		}

		nread += reader_npar(r, 0);
		reader_destroy(r);
	}
	idhash_build(psys.id, nread);
}

static void psystem_match_epoch(ReaderConstants * c) {
	/* match / read */
	int fid;
	intptr_t ipar;
	#pragma omp parallel for private(ipar)
	for(ipar = 0; ipar < psys.npar; ipar++) {
		psys.flag[ipar] |= PF_INVALID;
	}
	double distsum = 0.0;
	size_t lost = 0;
	for(fid = 0; fid < c->Nfiles; fid ++) {
		Reader * r = reader_new(psys.epoch->format);
		char * fname = reader_make_filename(psys.epoch->snapshot, fid);
		reader_open(r, fname);
		float (*pos)[3] = reader_alloc(r, "pos", 0);
		float * mass = reader_alloc(r, "mass", 0);
		float * sml = reader_alloc(r, "sml", 0);
		float * rho = reader_alloc(r, "rho", 0);
		float * ie = reader_alloc(r, "ie", 0);
		float * ye = reader_alloc(r, "ye", 0);
		float * xHI = reader_alloc(r, "xHI", 0); 
		/* FIXME: xHI not used yet, maybe a good idea to calculate
		 * the ye due to heavy elements and merge the contribution.
		 * */
		reader_read(r, "pos", 0, pos);
		reader_read(r, "mass", 0, mass);
		reader_read(r, "sml", 0, sml);
		reader_read(r, "rho", 0, rho);
		reader_read(r, "ie", 0, ie);
		reader_read(r, "ye", 0, ye);
		reader_read(r, "xHI", 0, xHI);

		uint64_t * id = calloc(8, reader_npar(r, 0));
		if(reader_itemsize(r, "id") == 4) {
			unsigned int * sid = reader_alloc(r, "id", 0);
			reader_read(r, "id", 0, id);
			intptr_t ipar;
			for(ipar = 0; ipar < reader_npar(r, 0); ipar++) {
				id[ipar] = sid[ipar];
			}
			free(sid);
		} else {
			reader_read(r, "id", 0, id);
		}

		intptr_t ipar;
		for(ipar = 0; ipar < reader_npar(r, 0); ipar++) {
			intptr_t hash = id[ipar] & IDHASHMASK;
			intptr_t best_jpar = -1;
			float mindist = 0.0;
			intptr_t jpar;
			for(jpar = psys.idhash.head[hash];
				jpar != -1; jpar = psys.idhash.next[jpar]) {
				float dst = dist(pos[ipar], psys.pos[jpar]);
				if(best_jpar == -1 || dst < mindist) {
					mindist = dst;
					best_jpar = jpar;
				}
			}
			if(best_jpar != -1) {
				memcpy(&psys.pos[best_jpar][0], &pos[ipar][0], sizeof(float) * 3);
				distsum += mindist;
				psys.mass[best_jpar] = mass[ipar];
				psys.sml[best_jpar] = sml[ipar];
				psys.rho[best_jpar] = rho[ipar];
				psys.ie[best_jpar] = ie[ipar];
				/*FIXME: somehow make use of the ye of the new snapshot*/
				psys.flag[best_jpar] &= ~PF_INVALID;
			} else {
				lost++;
			}
		}
		free(fname);
		free(mass);
		free(xHI);
		free(sml);
		free(rho);
		free(pos);
		free(ye);
		free(ie);
		free(id);
		reader_destroy(r);
	}
	MESSAGE("EPOCH matching: %ld skipped;  mean pos shifting = %lg", lost, 
		distsum/ c->Ntot[0]);
}

void solve_u_v(double d[3], double u[3], double v[3]) {
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

static void psystem_read_source() {
	FILE * f = fopen(psys.epoch->source, "r");
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

static float dist(const float p1[3], const float p2[3]) {
	int d;
	double result = 0.0;
	for (d = 0; d < 3; d++) {
		float dd = fabs(p1[d] - p2[d]);
		if(psys.periodic) {
			if (dd > 0.5 * psys.boxsize) dd = psys.boxsize - dd;
		}
		result += dd *dd;
	}
	return sqrt(result);
}

void psystem_write_output(int outputnum) {
	char * basename;
	asprintf(&basename, psys.epoch->output.filename, outputnum);
	int fid = 0;
	size_t nfiles = psys.epoch->output.nfiles;
	for(fid = 0; fid < nfiles; fid++) {
		char * filename;
		if(nfiles == 1) {
			filename = strdup(basename);
		} else {
			asprintf(&filename, "%s.%d", basename, fid);
		}
		/* write the gas */
		intptr_t gas_start = psys.npar * fid / nfiles;
		intptr_t gas_end = psys.npar * (fid + 1)/ nfiles;
		intptr_t gas_size = gas_end - gas_start;
		intptr_t bh_start = psys.nsrcs * fid / nfiles;
		intptr_t bh_end = psys.nsrcs * (fid + 1)/ nfiles;
		intptr_t bh_size = bh_end - bh_start;
		intptr_t i;
		Reader * r = reader_new("psphray");
		reader_create(r, filename);
		ReaderConstants * c = reader_constants(r);
		c->hasHe = 1;
		c->Ntot[0] = psys.npar;
		c->N[0] = gas_size;
		c->Ntot[5] = psys.nsrcs;
		c->N[5] = bh_size;

		c->Nfiles = nfiles;
		c->OmegaL = C_OMEGA_L;
		c->OmegaM = C_OMEGA_M;
		c->OmegaB = C_OMEGA_B;
		double znow = t2z(psys.epoch->age + psys.tick * psys.tick_time);
		c->time = 1. / (1. + znow);
		c->h = C_H;
		c->redshift = psys.epoch->redshift;
		c->boxsize = psys.boxsize;
		reader_update_header(r);
		reader_write(r, "pos", 0, psys.pos[gas_start]);
		reader_write(r, "sml", 0, &psys.sml[gas_start]);
		reader_write(r, "rho", 0, &psys.rho[gas_start]);
		reader_write(r, "mass", 0, &psys.mass[gas_start]);
		float * xHI = reader_alloc(r, "xHI", 0);
		float * ye = reader_alloc(r, "ye", 0);
		float * xHeI = reader_alloc(r, "xHeI", 0);
		float * xHeII = reader_alloc(r, "xHeII", 0);
		for(i = 0; i < gas_size; i++) {
			xHI[i] = psys_xHI(i);
			xHeI[i] = psys_xHeI(i);
			xHeII[i] = psys_xHeII(i);
			ye[i] = psys_ye(i);
		}
		reader_write(r, "xHI", 0, xHI);
		reader_write(r, "xHeI", 0, xHeI);
		reader_write(r, "xHeII", 0, xHeII);
		reader_write(r, "ye", 0, ye);
		free(xHeI);
		free(xHeII);
		free(xHI);
		free(ye);
		reader_write(r, "ie", 0, &psys.ie[gas_start]);
		reader_write(r, "lasthit", 0, &psys.lasthit[gas_start]);
		reader_write(r, "hits", 0, &psys.hits[gas_start]);
		reader_write(r, "id", 0, &psys.id[gas_start]);
		/* write the bh */
		float (*pos)[3] = reader_alloc(r, "pos", 5);
		for(i = 0; i < bh_size; i++) {
			pos[i][0] = psys.srcs[i+bh_start].pos[0];
			pos[i][1] = psys.srcs[i+bh_start].pos[1];
			pos[i][2] = psys.srcs[i+bh_start].pos[2];
		}
		reader_write(r, "pos", 5, pos);
		free(pos);
		double * ngammas = reader_alloc(r, "ngammas", 5);
		for(i = 0; i < bh_size; i++) {
			ngammas[i] = psys.srcs[i+bh_start].Ngamma_dot;
		}
		reader_write(r, "ngammas", 5, ngammas);
		free(ngammas);

		reader_close(r);
		reader_destroy(r);
		free(filename);
	}
	free(basename);
}

static inline void _get_value(void * field, intptr_t ipar, int type, int dim, double * value) {
	int d;
	for(d = 0; d < dim; d++) {
		switch(type) {
			case 0:
				value[d] = ((float*)field)[ipar * dim + d];
			break;
			case 1:
				value[d] = ((double*)field)[ipar * dim + d];
			break;
		}
	}
}

static void psystem_stat_internal(void * field, size_t npar, int type, int dim, double * max, double * min, double * mean) {
	intptr_t ipar = 0;
	double mn[dim], mx[dim], sum[dim];
	double value[dim];
	int d;
	_get_value(field, 0, type, dim, value);
	for(d = 0; d < dim; d++) {
		mn[d] = value[d]; 
		mx[d] = value[d];
		sum[d] = value[d];
	}

	#pragma omp parallel private(ipar)
	{
	double _mn[dim], _mx[dim], _sum[dim];
	double value[dim];
	int d;
	for(d = 0; d < dim; d++) {
		_mn[d] = mn[d]; 
		_mx[d] = mx[d];
		_sum[d] = 0.0;
	}
	#pragma omp for
	for(ipar = 1; ipar < npar; ipar++) {
		_get_value(field, ipar, type, dim, value);
		for(d = 0; d < dim; d++) {
			if(isnan(value[d])) {
				WARNING("nan found ipar=%ld", ipar);
			}
			if(_mx[d] < value[d]) _mx[d] = value[d];
			if(_mn[d] > value[d]) _mn[d] = value[d];
			_sum[d] += value[d];
		}
	}
	#pragma omp critical
	for(d = 0; d < dim; d++) {
		if(mx[d] < _mx[d]) mx[d] = _mx[d];
		if(mn[d] > _mn[d]) mn[d] = _mn[d];
		sum[d] += _sum[d];
	}
	}
	for(d = 0; d < dim; d++) {
		max[d] = mx[d];
		min[d] = mn[d];
		mean[d] = sum[d] / npar;
	}
}

void psystem_stat(const char * component) {
	double min[3], max[3], mean[3];
	if(!strcmp(component, "lambdaH")) {
		psystem_stat_internal(psys.lambdaH, psys.npar, 1, 1, max, min, mean);
	}
	if(!strcmp(component, "xHI")) {
		float * xHI = malloc(sizeof(float) * psys.npar);
		intptr_t i;
		#pragma omp parallel for private(i)
		for(i = 0; i < psys.npar; i++) {
			xHI[i] = psys_xHI(i);
		}
		psystem_stat_internal(xHI, psys.npar, 0, 1, max, min, mean);
		free(xHI);
	}
	if(!strcmp(component, "xHeI")) {
		float * xHeI = malloc(sizeof(float) * psys.npar);
		intptr_t i;
		#pragma omp parallel for private(i)
		for(i = 0; i < psys.npar; i++) {
			xHeI[i] = psys_xHeI(i);
		}
		psystem_stat_internal(xHeI, psys.npar, 0, 1, max, min, mean);
		free(xHeI);
	}
	if(!strcmp(component, "yGrecHII")) {
		psystem_stat_internal(psys.yGrecHII, psys.npar, 0, 1, max, min, mean);
	}
	if(!strcmp(component, "yGrecHeII")) {
		psystem_stat_internal(psys.yGrecHeII, psys.npar, 0, 1, max, min, mean);
	}
	if(!strcmp(component, "yGrecHeIII")) {
		psystem_stat_internal(psys.yGrecHeIII, psys.npar, 0, 1, max, min, mean);
	}
	if(!strcmp(component, "yGdepHI")) {
		psystem_stat_internal(psys.yGdepHI, psys.npar, 0, 1, max, min, mean);
	}
	if(!strcmp(component, "yGdepHeI")) {
		psystem_stat_internal(psys.yGdepHeI, psys.npar, 0, 1, max, min, mean);
	}
	if(!strcmp(component, "yGdepHeII")) {
		psystem_stat_internal(psys.yGdepHeII, psys.npar, 0, 1, max, min, mean);
	}
	if(!strcmp(component, "heat")) {
		psystem_stat_internal(psys.heat, psys.npar, 0, 1, max, min, mean);
	}
	if(!strcmp(component, "ye")) {
		float * ye = malloc(sizeof(float) * psys.npar);
		intptr_t i;
		#pragma omp parallel for private(i)
		for(i = 0; i < psys.npar; i++) {
			ye[i] = psys_ye(i);
		}

		psystem_stat_internal(ye, psys.npar, 0, 1, max, min, mean);
		free(ye);
	}
	if(!strcmp(component, "yeMET")) {
		psystem_stat_internal(psys.yeMET, psys.npar, 0, 1, max, min, mean);
	}
	if(!strcmp(component, "mass")) {
		psystem_stat_internal(psys.mass, psys.npar, 0, 1, max, min, mean);
	}
	if(!strcmp(component, "sml")) {
		psystem_stat_internal(psys.sml, psys.npar, 0, 1, max, min, mean);
	}
	if(!strcmp(component, "rho")) {
		psystem_stat_internal(psys.rho, psys.npar, 0, 1, max, min, mean);
	}
	if(!strcmp(component, "ie")) {
		psystem_stat_internal(psys.ie, psys.npar, 0, 1, max, min, mean);
	}
	if(!strcmp(component, "T")) {
		float * T = malloc(sizeof(float) * psys.npar);
		intptr_t i;
		#pragma omp parallel for private(i)
		for(i = 0; i < psys.npar; i++) {
			T[i] = psys_T(i);
		}
		psystem_stat_internal(T, psys.npar, 0, 1, max, min, mean);
		free(T);
	}
	MESSAGE("%s: mean=%g min=%g max=%g", component, mean[0], min[0], max[0]);
}
