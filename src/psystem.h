typedef struct {
	float pos[3];
	double Ngamma_sec;
	intptr_t lastemit;
	int specid;
} Source;

typedef struct {
	float (*pos)[3];
	double *xHI;
	double *ye;
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
