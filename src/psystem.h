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
	float *T;
	char * mask;
	size_t npar;
	int lasthit;
	float * deposit;
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
} PSystem;
