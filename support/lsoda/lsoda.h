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
	int jt;
	int mu;
	int ml;
};
typedef void    (*_lsoda_f) (double, double *, double *, void *);

void lsoda(_lsoda_f f, int neq, double *y, double *t, double tout, int itol, double *rtol, double *atol,
		   int itask, int *istate, struct lsoda_opt_t * opt, void *_data);
