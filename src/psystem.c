#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <messages.h>
#include <bitmask.h>

#include "config.h"
#include "reader.h"
#include "psystem.h"

#define IDHASHBITS 24
#define IDHASHMASK ((((size_t)1) << IDHASHBITS) - 1)

PSystem psys = {0};

void idhash_build(unsigned long long * id, size_t n) {
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

static float dist(const float p1[3], const float p2[3]);
void psystem_switch_epoch(int i) {


	Reader * r0 = reader_new(EPOCHS[i].format);
	char * fname0 = reader_make_filename(EPOCHS[i].snapshot, 0);
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
	MESSAGE("epoch %d, redshift %g, snapshot %s, format %s, source %s, ticks %u, ngas %lu, age %f [myr], duration %f [myr]",
	i, 
	EPOCHS[i].redshift,
	EPOCHS[i].snapshot,
	EPOCHS[i].format,
	EPOCHS[i].source,
	EPOCHS[i].nticks,
	EPOCHS[i].ngas,
	units_format(EPOCHS[i].age, "myr"),
	units_format(EPOCHS[i].duration, "myr"));

	if(i == 0) {
		size_t ngas = EPOCHS[i].ngas;
		psys.npar = ngas;
		int fid;
		size_t nread = 0;
		psys.pos = calloc(sizeof(float) * 3, ngas);
		psys.mass = calloc(sizeof(float), ngas);
		psys.id = calloc(sizeof(unsigned long long), ngas);
		psys.ie = calloc(sizeof(float), ngas);
		psys.sml = calloc(sizeof(float), ngas);
		psys.rho = calloc(sizeof(float), ngas);
		psys.xHI = calloc(sizeof(float), ngas);
		psys.ye = calloc(sizeof(float), ngas);
		psys.recomb = calloc(sizeof(float), ngas);
		psys.deposit = calloc(sizeof(float), ngas);
		psys.lasthit = calloc(sizeof(intptr_t), ngas);

		psys.mask = bitmask_alloc(ngas);

		for(fid = 0; fid < c->Nfiles; fid ++) {
			Reader * r = reader_new(EPOCHS[i].format);
			char * fname = reader_make_filename(EPOCHS[i].snapshot, fid);
			reader_open(r, fname);
			free(fname);
			reader_read(r, "pos", 0, psys.pos[nread]);
			reader_read(r, "mass", 0, &psys.mass[nread]);
			reader_read(r, "sml", 0, &psys.sml[nread]);
			reader_read(r, "rho", 0, &psys.rho[nread]);
			reader_read(r, "ye", 0, &psys.ye[nread]);
			reader_read(r, "ie", 0, &psys.ie[nread]);

			if(reader_itemsize(r, "id") == 4) {
				unsigned int * id = reader_alloc(r, "id", 0);
				reader_read(r, "id", 0, id);
				intptr_t ipar;
				for(ipar = 0; ipar < reader_npar(r, 0); ipar++) {
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
		bitmask_set_all(psys.mask);
	} else {
		int fid;
		bitmask_clear_all(psys.mask);
		double distsum = 0.0;
		size_t lost = 0;
		for(fid = 0; fid < c->Nfiles; fid ++) {
			Reader * r = reader_new(EPOCHS[i].format);
			char * fname = reader_make_filename(EPOCHS[i].snapshot, fid);
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

			unsigned long long * id = calloc(8, reader_npar(r, 0));
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
					bitmask_set(psys.mask, best_jpar);
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
		MESSAGE("EPOCH matching: %ld skipped;  mean pos shifting = %lg\n", lost, 
			distsum/ c->Ntot[0]);
	}
	MESSAGE("EPOCH active gas particles %ld/%ld\n", bitmask_sum(psys.mask), c->Ntot[0]);
	reader_destroy(r0);

	psys.tick = 0;
	psys.nticks = EPOCHS[i].nticks;
	psys.tick_time = EPOCHS[i].duration / EPOCHS[i].nticks;

	if(EPOCHS[i].source) {
		free(psys.srcs);
		FILE * f = fopen(EPOCHS[i].source, "r");
		if(f == NULL) {
			ERROR("failed to open %s", EPOCHS[i].source);
		}
		int NR = 0;
		char * line = NULL;
		size_t len = 0;
		intptr_t isrc = 0;
		int stage = 0;
		double x, y, z, dummy, L;
		char type[128];
		char garbage[128];

		while(0 <= getline(&line, &len, f)) {
			if(line[0] == '#') {
				NR++;
				continue;
			}
			switch(stage) {
			case 0:
				if(1 != sscanf(line, "%ld", &psys.nsrcs)) {
					ERROR("%s format error at %d", EPOCHS[i].source, NR);
				} else {
					psys.srcs = calloc(sizeof(Source), psys.nsrcs);
					isrc = 0;
					stage ++;
				}
				break;
			case 1:
				if(9 != sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %80s %80s",
					&x, &y, &z, &dummy, &dummy, &dummy, 
					&L, type, garbage)) {
					ERROR("%s format error at %d", EPOCHS[i].source, NR);
				}
				psys.srcs[isrc].pos[0] = x;
				psys.srcs[isrc].pos[1] = y;
				psys.srcs[isrc].pos[2] = z;
				psys.srcs[isrc].Ngamma_sec = L * 1e50;
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
			ERROR("%s expecting %lu, but has %lu sources, %lu", EPOCHS[i].source, psys.nsrcs, isrc, NR);
		}
		free(line);
		fclose(f);
	}

	intptr_t isrc;
	double Lmin;
	double Lmax;
	intptr_t imin = -1, imax = -1;
	for(isrc = 0; isrc < psys.nsrcs; isrc++) {
		if(isrc == 0 || psys.srcs[isrc].Ngamma_sec > Lmax) {
			Lmax = psys.srcs[isrc].Ngamma_sec;
			imax = isrc;
		}
		if(isrc == 0 || psys.srcs[isrc].Ngamma_sec < Lmin) {
			Lmin = psys.srcs[isrc].Ngamma_sec;
			imin = isrc;
		}
	}
	MESSAGE("SOURCES: %lu, max=%g at (%g %g %g), min=%g at (%g %g %g)",
		psys.nsrcs, Lmax, 
		psys.srcs[imax].pos[0],
		psys.srcs[imax].pos[1],
		psys.srcs[imax].pos[2],
		Lmin,
		psys.srcs[imin].pos[0],
		psys.srcs[imin].pos[1],
		psys.srcs[imin].pos[2]);
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

