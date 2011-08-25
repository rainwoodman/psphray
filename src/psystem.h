#include <math.h>
#ifndef MAXFLOAT
#include <float.h>
#define MAXFLOAT FLT_MAX
#endif
typedef struct {
	float pos[3];
	double Ngamma_sec;
	intptr_t lastemit;
	int specid;
} Source;

#define PF_NORMAL (0)
#define PF_INVALID (1)

typedef struct {
	float (*pos)[3];
	double *lambdaHI;  /* lambdaHI = arctan(xHI / xHII) */
	double *yeMET;
	float *mass;
	float *sml;
	float *rho;
	float *ie;
	int8_t * flag;
	size_t npar;

	double * recomb; /* number of photons from recombination */
	intptr_t * lasthit; /* time tick of last update */

	uint64_t * id;
	struct {
		intptr_t *head;
		intptr_t *next;
	} idhash;
	double boxsize;
	int periodic;
	Source * srcs;
	size_t nsrcs;
	intptr_t tick;
	double tick_time;
	Epoch * epoch;
} PSystem;

void psys_switch_epoch(int epoch);
void psystem_write_output(int outputnum);
void psystem_stat(const char * component);

typedef struct _Step {
	double lambdaHI;
	double xHI;
	double xHII;
	double ye;
	double ie;
	double yeMET;

	double nH;
	double T;

	double dxHI;
	double dye;
	double die;
	double dyGH;
} Step;

typedef struct x_t {
	intptr_t ipar;
	float d;
	float b;
} Xtype;

extern PSystem psys;

static inline const double psys_xHI(const intptr_t i) {
	const double r = tan(psys.lambdaHI[i]);
	return r / (1. + r);
}

#define lambdaHI_to_xHI_xHII(lambda, xHI, xHII) \
{ const double r = tan(lambda); \
  const double r_inv = (r!=0.0)?1.0/r:MAXFLOAT ; \
  xHI = r / ( 1. + r); \
  xHII = r_inv / (1. + r_inv); \
  if(fabs(xHI + xHII - 1.0) > 1e-1) abort(); \
}

static inline const double psys_xHII(const intptr_t i) {
	const double r = tan(psys.lambdaHI[i]);
	const double r_inv = (r!=0.0)?1/r:MAXFLOAT;
	return r_inv / (1. + r_inv);
}

static inline const double psys_ye(const intptr_t i) {
	return psys_xHII(i) + psys.yeMET[i];
}
static inline const double psys_T(const intptr_t i) {
	return ieye2T(psys.ie[i], psys_ye(i));
}
static inline void psys_set_lambdaHI(const intptr_t i, const double xHI, const double xHII) {
	if(xHI >= 1.0 || xHII <= 0.0) psys.lambdaHI[i] = atan(MAXFLOAT);
	else if(xHI <= 0.0 || xHII >=1.0) psys.lambdaHI[i] = 0.0;
	else psys.lambdaHI[i] = atan2(xHI , xHII);
}

#define __psystem_swap__(a, b, type) {const type __tmp__ = a; a = b; b = __tmp__;}
/* swapping two particles note that IDHash will be corrupted after this call */
static inline void psys_swap(intptr_t i, intptr_t j) {
	int d;
	for (d = 0; d < 3; d++) {
		__psystem_swap__(psys.pos[i][d], psys.pos[j][d], float);
	}
	__psystem_swap__(psys.lambdaHI[i], psys.lambdaHI[j], double);
	__psystem_swap__(psys.yeMET[i], psys.yeMET[j], double);
	__psystem_swap__(psys.mass[i], psys.mass[j], float);
	__psystem_swap__(psys.sml[i], psys.sml[j], float);
	__psystem_swap__(psys.rho[i], psys.rho[j], float);
	__psystem_swap__(psys.ie[i], psys.ie[j], float);
	__psystem_swap__(psys.flag[i], psys.flag[j], uint8_t);
	__psystem_swap__(psys.recomb[i], psys.recomb[j], double);
	__psystem_swap__(psys.lasthit[i], psys.lasthit[j], intptr_t);
	__psystem_swap__(psys.id[i], psys.id[j], uint64_t);
}
