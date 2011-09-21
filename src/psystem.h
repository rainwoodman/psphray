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

	double ray_length_hint;

	/* below for type = SRC_PLANE only*/
	double radius;
	double a[3];
	double b[3];
} Source;

#define PF_NORMAL (0)
#define PF_INVALID (1)
#define PF_HOTSPOT (8)

typedef struct {
	float (*pos)[3];
	double *lambdaH;  /* lambdaH = arctan(xHI / xHII) */
	double *lambdaHeI;  /* lambdaHeI = arctan(xHeI / (xHeII + xHeIII)) */
	double *lambdaHeII;  /* lambdaHeII = arctan(xHeII / xHeIII) */
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
	int * hits;   /* total number of hits, subtotaled between snapshots*/

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
	double lambdaH;
	double ie;
	double yeMET;
	double yGdep;
	double heat;
	double time;
	double nH;
	double T;

	double dyGrec;
} Step;


typedef struct x_t {
	intptr_t ipar;
	double d;
	double b;
} Xtype;

extern PSystem psys;

/* conversion from lambdaH to xHI and xHII */
static inline const double lambdaH_to_xHI(const double lambdaH) {
	const double r = tan(lambdaH);
	return r / (1. + r);
}
static inline const double lambdaH_to_xHII(const double lambdaH) {
	const double r = tan(lambdaH);
	const double r_inv = (r!=0.0)?1/r:HUGENUMBER;
	return r_inv / (1. + r_inv);
}
static inline const double lambdaH_from_xHI_xHII(const double xHI, const double xHII) {
}
/* conversion from lambdaHeI, lambdaHeII to xHeI and xHeII and xHeIII */
static inline const double lambdaHe_to_xHeI(const double lambdaHeI, const double lambdaHeII) {
	const double r = tan(lambdaHeI);
	return r / (1. + r);
}
static inline const double lambdaHe_to_xHeII(const double lambdaHeI, const double lambdaHeII) {
	const double r1 = tan(lambdaHeI);
	const double r2 = tan(lambdaHeII);
	const double r1_inv = (r1!=0.0)?1/r1:HUGENUMBER;
	return r1_inv / (1. + r1_inv) * r2 / (1. + r2);
}

static inline const double lambdaHe_to_xHeIII(const double lambdaHeI, const double lambdaHeII) {
	const double r1 = tan(lambdaHeI);
	const double r2 = tan(lambdaHeII);
	const double r1_inv = (r1!=0.0)?1/r1:HUGENUMBER;
	const double r2_inv = (r2!=0.0)?1/r2:HUGENUMBER;
	return r1_inv / (1. + r1_inv) * r2_inv / (1. + r2_inv);
}


static inline const double psys_xHI(const intptr_t i) {
	return lambdaH_to_xHI(psys.lambdaH[i]);
}
static inline const double psys_xHII(const intptr_t i) {
	return lambdaH_to_xHII(psys.lambdaH[i]);
}
static inline const double psys_xHeI(const intptr_t i) {
	return lambdaHe_to_xHeI(psys.lambdaHeI[i], psys.lambdaHeII[i]);
}
static inline const double psys_xHeII(const intptr_t i) {
	return lambdaHe_to_xHeII(psys.lambdaHeI[i], psys.lambdaHeII[i]);
}
static inline const double psys_xHeIII(const intptr_t i) {
	return lambdaHe_to_xHeIII(psys.lambdaHeI[i], psys.lambdaHeII[i]);
}

static inline void psys_set_lambdaH(const intptr_t i, const double xHI, const double xHII) {
	if(xHI >= 1.0 || xHII <= 0.0) psys.lambdaH[i] = atan(HUGENUMBER);
	else if(xHI <= 0.0 || xHII >=1.0) psys.lambdaH[i] = atan(1/HUGENUMBER);
	else psys.lambdaH[i] = atan2(xHI , xHII);
}
static inline void psys_set_lambdaHe(const intptr_t i, const double xHeI, const double xHeII, const double xHeIII) {
	if(xHeI >= 1.0 || (xHeII + xHeIII) <= 0.0) psys.lambdaHeI[i] = atan(HUGENUMBER);
	else if(xHeI <= 0.0 || (xHeII + xHeIII) >=1.0) psys.lambdaHeI[i] = atan(1/HUGENUMBER);
	else psys.lambdaHeI[i] = atan2(xHeI , xHeII + xHeIII);

	if((xHeII + xHeI) >= 1.0 || xHeIII <= 0.0) psys.lambdaHeII[i] = atan(HUGENUMBER);
	else if(xHeII <= 0.0 || (xHeIII + xHeI) >=1.0) psys.lambdaHeII[i] = atan(1/HUGENUMBER);
	else psys.lambdaHeII[i] = atan2(xHeII , xHeIII);
}

static inline const double psys_NH(const intptr_t i) {
	return psys.mass[i] * C_H_PER_MASS;
}
static inline const double psys_NHe(const intptr_t i) {
	return psys.mass[i] * C_HE_PER_MASS;
}
static inline const double psys_nH(const intptr_t i) {
	return psys.rho[i] * C_H_PER_MASS;
}
static inline const double psys_nHe(const intptr_t i) {
	return psys.rho[i] * C_HE_PER_MASS;
}

static inline const double psys_ye(const intptr_t i) {
	return psys_xHII(i) + psys_xHeII(i) + 2. * psys_xHeIII(i) + psys.yeMET[i];
}
static inline const double psys_T(const intptr_t i) {
	return ieye2T(psys.ie[i], psys_ye(i));
}
