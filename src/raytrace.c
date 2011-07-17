#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <messages.h>
#include <bitmask.h>
#include "psystem.h"

typedef struct _Cell Cell;
struct _Cell {
	intptr_t head_par;
	intptr_t first_child;
	intptr_t parent;
	float bot[3];
	float top[3];
	unsigned int npar;
};

typedef struct {
	intptr_t * next;
	Cell * pool;
	intptr_t pool_used;
	size_t pool_size;
} Raytrace;

extern PSystem psys;

extern int pluecker_(float dir[3], float dist, float s2b[3], float s2t[3]);
static int inside(float pos[3], intptr_t icell);
static int full(intptr_t icell);
static intptr_t find(float pos[3], intptr_t icell);
static intptr_t parent(intptr_t icell);
static void add(intptr_t ipar, intptr_t icell);
static int split(intptr_t icell);

Raytrace rt = {0};
const size_t PPC = 16;

void rt_switch_epoch(int i) {
	if(i == 0) {
		rt.next = malloc(sizeof(intptr_t) * psys.npar);
		rt.pool_size = psys.npar / PPC;
		rt.pool = malloc(sizeof(Cell) * rt.pool_size);
	}
tryagain:
	rt.pool_used = 0;
	memset(rt.next, -1, sizeof(intptr_t) * psys.npar);
	rt.pool[0].parent = -1;
	int d;
	for(d = 0; d < 3; d++) {
		rt.pool[0].bot[d] = 0;
		rt.pool[0].top[d] = psys.boxsize;
	}
	rt.pool[0].first_child = -1;
	rt.pool[0].head_par = -1;
	rt.pool[0].npar = 0;
	rt.pool_used ++;

	intptr_t icell = 0;
	intptr_t ipar;
	size_t skipped = 0;
	for(ipar = 0; ipar < psys.npar; ipar++) {
		float * pos = psys.pos[ipar];
		while(!inside(pos, icell)) {
			icell = parent(icell);
			if(icell == -1) break;
		}
		/* if par not in the root cell, skip it*/
		if(icell == -1) {
			skipped ++;
			continue;
		}

		icell = find(pos, icell);

		while(full(icell)) {
			if(!split(icell)) {
				MESSAGE("OCTTREE full %lu mean occupy = %f", rt.pool_size, (double)ipar / rt.pool_size);
				size_t trysize = psys.npar / ((double) ipar / rt.pool_size);
				free(rt.pool);
				if(rt.pool_size < trysize) rt.pool_size = trysize;
				else rt.pool_size *= 2;
				rt.pool = malloc(sizeof(Cell) * rt.pool_size);
				MESSAGE("retry with %lu", rt.pool_size);
				goto tryagain;
			}
			icell = find(pos, icell);
		}
		add(ipar, icell);
	}
	if(skipped > 0)
		WARNING("%ld particles out of cell", skipped);
	MESSAGE("Cells %lu/%lu, mean occupicy = %f", rt.pool_used, rt.pool_size, 
			(double)psys.npar / rt.pool_used);

	/* now recalculate the AABB boxes */
	size_t done = 0;
	char * child_done = malloc(rt.pool_used);
	memset(child_done, 0, rt.pool_used);

	for(icell = 0; icell < rt.pool_used; icell++) {
		if(rt.pool[icell].first_child == -1) {
			child_done[icell] = 8;
		}
	}
	
	while(done < rt.pool_used) {	
		for(icell = 0; icell < rt.pool_used; icell++) {
			if(child_done[icell] > 8) ERROR("more than 8 children?");
			if(child_done[icell] != 8) continue;
			float * top = rt.pool[icell].top;
			float * bot = rt.pool[icell].bot;
			if(rt.pool[icell].first_child != -1) {
				int first = 1;
				for(ipar = rt.pool[icell].head_par; ipar!= -1; ipar = rt.next[ipar]) {
					float * pos = psys.pos[ipar];
					float sml = psys.sml[ipar];
					for(d = 0; d < 3; d++) {
						if(first || pos[d] - sml < bot[d]) bot[d] = pos[d] - sml;
						if(first || pos[d] + sml > top[d]) top[d] = pos[d] + sml;
					}
					first = 0;
				}
			} else {
				int first = 1;
				int i;
				for(i = 0; i < 8 ;i++) {
					intptr_t first_child = rt.pool[icell].first_child;
					float * cbot = rt.pool[first_child + i].bot;
					float * ctop = rt.pool[first_child + i].top;
					for(d = 0; d < 3; d++) {
						if(first || cbot[d] < bot[d]) bot[d] = cbot[d];
						if(first || ctop[d] > top[d]) top[d] = ctop[d];
					}
					first = 0;
				}
			}
			done ++;
			child_done[icell] = 0;
			if(icell != 0) {
				child_done[rt.pool[icell].parent]++;
			}
		}
		MESSAGE("updating AABB %lu/%lu done ", done, rt.pool_used);
	}
	free(child_done);
}

static int hit(float s[3], float dir[3], float dist, intptr_t icell) {
	float * bot = rt.pool[icell].bot;
	float * top = rt.pool[icell].top;
	float s2b[3], s2t[3];
	int d;
	for(d = 0; d < 3; d++) {
		s2b[d] = bot[d] - s[d];
		s2t[d] = top[d] - s[d];
	}
	return pluecker_(dir, dist, s2b, s2t);
}

static void add(intptr_t ipar, intptr_t icell) {
	Cell * cell = &rt.pool[icell];
	if(cell->first_child != -1) {
		ERROR("never shall reach here, adding to a none leaf")
	}
	rt.next[ipar] = cell->head_par;
	cell->head_par = ipar;
	cell->npar++;
}
static int full(intptr_t icell) {
	return rt.pool[icell].npar >= PPC;
}
static intptr_t parent(intptr_t icell) {
	return rt.pool[icell].parent;
}
static intptr_t find(float pos[3], intptr_t icell) {
	Cell * cell = &rt.pool[icell];
	if(cell->first_child == -1) return icell;
	int i;
	for(i = 0; i < 8; i++) {
		if(inside(pos, cell->first_child+i)) {
			return find(pos, cell->first_child + i);
		}
	}
	ERROR("never reach here makesure inside() is ensured before find()");
}
static int split(intptr_t icell) {
	static int bitmask[] = {0x1, 0x2, 0x4};
	Cell * cell = &rt.pool[icell];
	if(rt.pool_used +8 >=rt.pool_size) {
		return 0;
	}
	cell->first_child = rt.pool_used;
	rt.pool_used += 8;
	Cell * fc = &rt.pool[cell->first_child];
	float center[3];
	int d;
	int i;
	for(d = 0; d < 3; d++) {
		center[d]  = 0.5 * (cell->top[d] + cell->bot[d]);
	}
	for(i = 0; i < 8; i++) {
		fc[i].parent = icell;
		for(d = 0; d < 3; d++) {
			fc[i].bot[d] = (bitmask[d] & i)?center[d]:cell->bot[d];
			fc[i].top[d] = (bitmask[d] & i)?cell->top[d]:center[d];
		}
		fc[i].head_par = -1;
		fc[i].first_child = -1;
		fc[i].npar = 0;
	}
	cell->npar = 0;
	intptr_t ipar;
	intptr_t nextpar;
	for(ipar = cell->head_par; ipar != -1; ipar = nextpar) {
		nextpar = rt.next[ipar];
		for(i = 0; i < 8; i++) {
			if(inside(psys.pos[ipar], i + cell->first_child)) {
				add(ipar, i + cell->first_child);
			}
		}
	}
	cell->head_par = -1;
	return 1;
}
static inside(float pos[3], intptr_t icell) {
	int d;
	Cell * cell = &rt.pool[icell];
	for(d = 0; d < 3; d++) {
		if(pos[d] >= cell->top[d] ||
			pos[d] < cell->bot[d]) return 0;
	}
	return 1;
}
