typedef struct {
	float pos[3];
	double Ngamma_sec;
} Source;

typedef struct {
	float (*pos)[3];
	float *xHI;
	float *ye;
	float *mass;
	float *sml;
	float *rho;
	float *ie;

	char * mask;
	size_t npar;

	float * recomb; /* number of photons from recombination */
	float * deposit; /* number of photons deposited since last update*/
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
	size_t nticks;
	double tick_time;
} PSystem;

void psys_switch_epoch(int epoch);
