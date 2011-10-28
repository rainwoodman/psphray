
#include <stdio.h>              /* used for NULL, also fprintf(stderr,...) */
#include <stdlib.h>             /* used for malloc's */
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>
#include "ompdiscrete.h"
#include "discrete.h"

ran_discrete_omp_t * ran_discrete_omp_preproc (size_t Kevents, const float *ProbArray) {
	ran_discrete_omp_t * g = calloc(1, sizeof(ran_discrete_omp_t));
	float * ptotal = NULL;

	#pragma omp parallel 
	{
		#pragma omp master
		{
			g->nthreads = omp_get_num_threads();
			g->rans = calloc(g->nthreads, sizeof(gsl_ran_discretef_t *));
			g->starts = calloc(g->nthreads, sizeof(size_t));
			ptotal = calloc(g->nthreads, sizeof(float));
		}
		#pragma omp barrier

		int tid = omp_get_thread_num();

		size_t start = Kevents * tid / g->nthreads;
		size_t end = Kevents * (tid + 1)/ g->nthreads;
		size_t len = end - start;
		g->starts[tid] = start;
		g->rans[tid] = gsl_ran_discretef_preproc0(len, &ProbArray[start], &ptotal[tid]);
		#pragma omp barrier
		#pragma omp master 
		{
			g->ranroot = gsl_ran_discretef_preproc(g->nthreads, ptotal);
			free(ptotal);
		}
	}
	return g;
}

size_t
ran_discrete_omp(const gsl_rng *r, const ran_discrete_omp_t *g) {
	size_t block = gsl_ran_discretef(r, g->ranroot);
	return g->starts[block] + gsl_ran_discretef(r, g->rans[block]);
}

void ran_discrete_omp_free(ran_discrete_omp_t * g) {
	if(g == NULL) return;
	int i;
	for(i = 0; i < g->nthreads; i++) {
		gsl_ran_discretef_free(g->rans[i]);
	}
	gsl_ran_discretef_free(g->ranroot);
	free(g->starts);
	free(g->rans);
	free(g);
}
