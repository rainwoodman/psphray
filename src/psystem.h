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
#define PSYS(field,i) (psys.field [i])

typedef struct {
	float (*pos)[3];
	double *lambdaH;  /* lambdaH = xHI) */
	double *lambdaHeI;  /* lambdaHeI = xHeI */
	double *lambdaHeII;  /* lambdaHeII = xHeII) */
	float *yeMET;
	float *mass;
	float *sml;
	float *rho;
	float *ie;
	int8_t * flag;
	size_t npar;
	size_t npar_max;

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

	double step_remain; /* how much (0. - 1.0) of the step remains after the solver */
	int refined; /* if the time step is too large and a refining is done in the intergrator */
	intptr_t ipar; /* debugging only the ipar used in this step */
} Step;


typedef struct x_t {
	intptr_t ipar;
	double d;
	double b;
} Xtype;

extern PSystem psys;

/* conversion from lambdaH to xHI and xHII */
static inline double lambdaH_to_xHI(const double lambdaH) {
	return lambdaH;
}
static inline double lambdaH_to_xHII(const double lambdaH) {
	return 1.0 - lambdaH;
}
/* conversion from lambdaHeI, lambdaHeII to xHeI and xHeII and xHeIII */
static inline double lambdaHe_to_xHeI(const double lambdaHeI, const double lambdaHeII) {
	return lambdaHeI;
}
static inline double lambdaHe_to_xHeII(const double lambdaHeI, const double lambdaHeII) {
	return lambdaHeII;
}

static inline double lambdaHe_to_xHeIII(const double lambdaHeI, const double lambdaHeII) {
	return 1.0 - lambdaHeII - lambdaHeI;
}
static inline double lambda_to_ye(const double lambdaH, double lambdaHeI, double lambdaHeII) {
	return lambdaH_to_xHII(lambdaH) 
		+ C_HEMF / C_HMF *0.25 * (lambdaHe_to_xHeII(lambdaHeI, lambdaHeII) + 2 * lambdaHe_to_xHeIII(lambdaHeI, lambdaHeII));

}

static inline double psys_xHI(const intptr_t i) {
	return lambdaH_to_xHI(PSYS(lambdaH, i));
}
static inline double psys_xHII(const intptr_t i) {
	return lambdaH_to_xHII(PSYS(lambdaH, i));
}
static inline double psys_xHeI(const intptr_t i) {
	return lambdaHe_to_xHeI(PSYS(lambdaHeI, i), PSYS(lambdaHeII, i));
}
static inline double psys_xHeII(const intptr_t i) {
	return lambdaHe_to_xHeII(PSYS(lambdaHeI, i), PSYS(lambdaHeII, i));
}
static inline double psys_xHeIII(const intptr_t i) {
	return lambdaHe_to_xHeIII(PSYS(lambdaHeI, i), PSYS(lambdaHeII, i));
}

static inline void psys_set_lambdaH(const intptr_t i, const double xHI, const double xHII) {
	if(xHI >= 1.0 || xHII <= 0.0) PSYS(lambdaH, i) = 1.0;
	else if(xHI <= 0.0 || xHII >=1.0) PSYS(lambdaH, i) = 0.0;
	else PSYS(lambdaH, i) = xHI;
}
static inline void psys_set_lambdaHe(const intptr_t i, const double xHeI, const double xHeII, const double xHeIII) {
	if(xHeI >= 1.0 || (xHeII + xHeIII) <= 0.0) PSYS(lambdaHeI, i) = 1.0;
	else if(xHeI <= 0.0 || (xHeII + xHeIII) >=1.0) PSYS(lambdaHeI, i) = 0.0;
	else PSYS(lambdaHeI, i) = xHeI;

	if(xHeII >= 1.0 || (xHeI + xHeIII) <= 0.0) PSYS(lambdaHeII, i) = 1.0;
	else if(xHeII <= 0.0 || (xHeIII + xHeI) >=1.0) PSYS(lambdaHeII, i) = 0.0;
	else PSYS(lambdaHeII, i) = xHeII;
}
static inline double psys_NH(const intptr_t i) {
	return PSYS(mass, i) * C_H_PER_MASS;
}
static inline double psys_NHe(const intptr_t i) {
	return PSYS(mass, i) * C_HE_PER_MASS;
}

static inline double psys_ye(const intptr_t i) {
	return lambda_to_ye(PSYS(lambdaH, i), PSYS(lambdaHeI, i), PSYS(lambdaHeII, i));
}
static inline double psys_T(const intptr_t i) {
	return ieye2T(PSYS(ie, i), psys_ye(i));
}

static inline double psys_Ngamma_dot(intptr_t i) {
	return PSYS(srcs, i).Ngamma_dots[PSYS(srcs, i).cursor];
}
void psystem_weight_srcs(double weights[]);

void psystem_get_source_weights(double weights[]);
