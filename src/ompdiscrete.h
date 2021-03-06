#include "discrete.h"

typedef struct {
	gsl_ran_discretef_t ** rans;
	gsl_ran_discretef_t * ranroot;
	int nthreads;
	size_t * starts;
} ran_discrete_omp_t;

ran_discrete_omp_t * ran_discrete_omp_preproc (size_t Kevents, const float *ProbArray);
size_t ran_discrete_omp(const gsl_rng *r, const ran_discrete_omp_t *g);
void ran_discrete_omp_free(ran_discrete_omp_t *g);


