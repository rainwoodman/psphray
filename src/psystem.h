#include <math.h>
#define HUGENUMBER 1e15
#define PSYS_SRC_PLANE (0)
#define PSYS_SRC_POINT (1)

typedef struct {
	double pos[3];
	double dir[3];
	double Ngamma_sec;
	intptr_t lastemit;
	int specid;
	int type;
	/* below for type = SRC_PLANE only*/
	double radius;
	double a[3];
	double b[3];
} Source;

#define PF_NORMAL (0)
#define PF_INVALID (1)

typedef struct {
	float (*pos)[3];
	double *lambdaHI;  /* lambdaHI = arctan(xHI / xHII) */
	float *yeMET;
	float *mass;
	float *sml;
	float *rho;
	float *ie;
	int8_t * flag;
	size_t npar;

	float * yGrec; /* number of recombination photon / NH */
	float * yGdep; /* number of ionizating photon / NH */
	float * heat;  /* total heat deposition per mass in a step*/
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
	double ie;
	double yeMET;
	double yGdep;
	double heat;
	double time;
	double nH;
	double rho;
	double T;

	double dyGrec;
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
  const double r_inv = (r!=0.0)?1.0/r:HUGENUMBER ; \
  xHI = r / ( 1. + r); \
  xHII = r_inv / (1. + r_inv); \
  if(fabs(xHI + xHII - 1.0) > 1e-1) abort(); \
}

static inline const double psys_NH(const intptr_t i) {
	return psys.mass[i] * C_H_PER_MASS;
}
static inline const double psys_nH(const intptr_t i) {
	return psys.rho[i] * C_H_PER_MASS;
}

static inline const double psys_xHII(const intptr_t i) {
	const double r = tan(psys.lambdaHI[i]);
	const double r_inv = (r!=0.0)?1/r:HUGENUMBER;
	return r_inv / (1. + r_inv);
}

static inline const double psys_ye(const intptr_t i) {
	return psys_xHII(i) + psys.yeMET[i];
}
static inline const double psys_T(const intptr_t i) {
	return ieye2T(psys.ie[i], psys_ye(i));
}
static inline void psys_set_lambdaHI(const intptr_t i, const double xHI, const double xHII) {
	if(xHI >= 1.0 || xHII <= 0.0) psys.lambdaHI[i] = atan(HUGENUMBER);
	else if(xHI <= 0.0 || xHII >=1.0) psys.lambdaHI[i] = atan(1/HUGENUMBER);
	else psys.lambdaHI[i] = atan2f(xHI , xHII);
}
