#include <math.h>
typedef struct {
	float pos[3];
	double Ngamma_sec;
	intptr_t lastemit;
	int specid;
} Source;

typedef struct {
	float (*pos)[3];
	double *lambdaHI;  /* lambdaHI = arctan(xHI / xHII) */
	double *yeMET;
	float *mass;
	float *sml;
	float *rho;
	float *ie;
	char * mask;
	size_t npar;

	double * recomb; /* number of photons from recombination */
	intptr_t * lasthit; /* time tick of last update */

	unsigned long long * id;
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
	double xHI;
	double xHII;
	double ye;
	double ie;
	double y;

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

static inline const double psys_xHII(const intptr_t i) {
	const double r = tan(psys.lambdaHI[i]);
	return 1.0 / (1. + r);
}

static inline const double psys_ye(const intptr_t i) {
	return psys_xHII(i) + psys.yeMET[i];
}
static inline const double psys_T(const intptr_t i) {
	return ieye2T(psys.ie[i], psys_ye(i));
}
static inline void psys_set_lambdaHI(const intptr_t i, const double xHI, const double xHII) {
	if(xHII > 1.0) psys.lambdaHI[i] = atan(HUGE_VALF);
	else if(xHI < 0.0) psys.lambdaHI[i] = 0.0;
	else psys.lambdaHI[i] = atan2(xHI, xHII);
}
