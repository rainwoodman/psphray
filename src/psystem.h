#include <math.h>
#include <array.h>
#define HUGENUMBER 1e15
#define PSYS_SRC_PLANE (0)
#define PSYS_SRC_POINT (1)

typedef struct {
	double pos[3];
	double dir[3];
	intptr_t lastemit;
	int specid;
	int type;

	double ray_length_hint;

	/* below for type = SRC_PLANE only*/
	double radius;
	double a[3];
	double b[3];

	/* luminosity history */
	ARRAY_DEFINE_S(ticks, intptr_t);
	ARRAY_DEFINE_S(Ngamma_dots, double);
	intptr_t cursor;
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

	union {
	struct {
	float * yGdepHI; /* number of ionizating photon / NH */
	float * yGdepHeI; /* number of ionizating photon / NHe */
	float * yGdepHeII; /* number of ionizating photon / NHe */
	};
	float * (yGdep[3]);
	};
	union {
	struct {
		float * yGrecHII; /* number of ionizating photon / NH */
		float * yGrecHeII; /* number of ionizating photon / NHe */
		float * yGrecHeIII; /* number of ionizating photon / NHe */
	};
	float * (yGrec[3]);
	};
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
	double lambdaHeI;
	double lambdaHeII;
	double ie;
	double yeMET;
	double yGdepHI;
	double yGdepHeI;
	double yGdepHeII;
	double heat;
	double time;
	double nH;
	double T;

	double yGrecHII;
	double yGrecHeII;
	double yGrecHeIII;
} Step;


typedef struct x_t {
	intptr_t ipar;
	double d;
	double b;
} Xtype;

extern PSystem psys;

/* conversion from lambdaH to xHI and xHII */
static inline const double lambdaH_to_xHI(const double lambdaH) {
	return lambdaH;
}
static inline const double lambdaH_to_xHII(const double lambdaH) {
	return 1.0 - lambdaH;
}
static inline const double lambdaH_from_xHI_xHII(const double xHI, const double xHII) {
}
/* conversion from lambdaHeI, lambdaHeII to xHeI and xHeII and xHeIII */
static inline const double lambdaHe_to_xHeI(const double lambdaHeI, const double lambdaHeII) {
	return lambdaHeI;
}
static inline const double lambdaHe_to_xHeII(const double lambdaHeI, const double lambdaHeII) {
	return 1.0 - lambdaHeI - lambdaHeII;
}

static inline const double lambdaHe_to_xHeIII(const double lambdaHeI, const double lambdaHeII) {
	return lambdaHeII;
}
static inline const double lambda_to_ye(const double lambdaH, double lambdaHeI, double lambdaHeII) {
	return lambdaH_to_xHII(lambdaH) 
		+ C_HEMF / C_HMF *0.25 * (lambdaHe_to_xHeII(lambdaHeI, lambdaHeII) + 2 * lambdaHe_to_xHeIII(lambdaHeI, lambdaHeII));

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
	if(xHI >= 1.0 || xHII <= 0.0) psys.lambdaH[i] = 1.0;
	else if(xHI <= 0.0 || xHII >=1.0) psys.lambdaH[i] = 0.0;
	else psys.lambdaH[i] = xHI;
}
static inline void psys_set_lambdaHe(const intptr_t i, const double xHeI, const double xHeII, const double xHeIII) {
	if(xHeI >= 1.0 || (xHeII + xHeIII) <= 0.0) psys.lambdaHeI[i] = 1.0;
	else if(xHeI <= 0.0 || (xHeII + xHeIII) >=1.0) psys.lambdaHeI[i] = 0.0;
	else psys.lambdaHeI[i] = xHeI;

	if(xHeIII >= 1.0 || (xHeI + xHeII) <= 0.0) psys.lambdaHeII[i] = 1.0;
	else if(xHeIII <= 0.0 || (xHeII + xHeI) >=1.0) psys.lambdaHeII[i] = 0.0;
	else psys.lambdaHeII[i] = xHeIII;
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
	return lambda_to_ye(psys.lambdaH[i], psys.lambdaHeI[i], psys.lambdaHeII[i]);
}
static inline const double psys_T(const intptr_t i) {
	return ieye2T(psys.ie[i], psys_ye(i));
}

static inline double psys_Ngamma_dot(intptr_t i) {
	return psys.srcs[i].Ngamma_dots[psys.srcs[i].cursor];
}
void psystem_weight_srcs(double weights[]);

void psystem_get_source_weights(double weights[]);
