int stoda(struct common_t * common, struct lsoda_context_t * ctx, double *y, int jstart, struct lsoda_opt_t * opt);
int      correction(struct common_t * common, struct lsoda_context_t * ctx, double *y, double pnorm, double *del, double *delp, double *told,
						   int *ncf, double *rh, int *m, double hmin);
int prja(struct common_t * common, struct lsoda_context_t * ctx, double *y);

int      ewset(struct common_t * common, const int neq, double ewt[], const double rtol[], const double atol[], const double *ycur);
void     solsy(struct common_t * common, int neq, double *y);
int      orderswitch(struct common_t * common, int neq, double rhup, double dsm, double *pdh, double *rh, int kflag, int maxord);
int intdy(struct common_t * common, int neq, double t, int k, double *dky);
int      corfailure(struct common_t * common, int neq, double *told, double *rh, int *ncf, double hmin);
void     methodswitch(struct common_t * common, int neq, double dsm, double pnorm, double *pdh, double *rh, int mxords, int mxordn);
void     cfode(struct common_t * common, int meth);
void     cfode_static(struct common_t * common, int meth);
void     scaleh(struct common_t * common, int neq, double *rh, double *pdh, double hmxi);

#ifdef CFODE_STATIC
	#define cfode cfode_static
#endif

