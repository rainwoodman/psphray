struct lsoda_opt_t {
	int ixpr;
	int mxstep;
	int mxhnil;
	int mxordn;
	int mxords;
	double tcrit;
	double h0;
	double hmax;
	double hmin;
	double hmxi;
};
typedef void    (*_lsoda_f) (double, double *, double *, void *);


struct lsoda_context_t {
	_lsoda_f function;
	void * data;
	int neq;
};


void lsoda(struct lsoda_context_t * ctx, double *y, double *t, double tout, int itol, double *rtol, double *atol,
		   int itask, int *istate, struct lsoda_opt_t * opt);
