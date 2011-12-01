int stoda(struct lsoda_context_t * ctx, double *y, int jstart);
int      correction(struct lsoda_context_t * ctx, double *y, double pnorm, double *del, double *delp, double *told,
						   int *ncf, double *rh, int *m, double hmin);
int prja(struct lsoda_context_t * ctx, double *y);

int      ewset(struct lsoda_context_t * ctx, const int neq, double ewt[], const double rtol[], const double atol[], const double *ycur);
void     solsy(struct lsoda_context_t * ctx, int neq, double *y);
int      orderswitch(struct lsoda_context_t * ctx, int neq, double rhup, double dsm, double *pdh, double *rh, int kflag, int maxord);
int      intdy(struct lsoda_context_t * ctx, int neq, double t, int k, double *dky);
int      corfailure(struct lsoda_context_t * ctx, int neq, double *told, double *rh, int *ncf, double hmin);
void     methodswitch(struct lsoda_context_t * ctx, int neq, double dsm, double pnorm, double *pdh, double *rh, int mxords, int mxordn);
void     cfode(struct lsoda_context_t * ctx, int meth);
void     cfode_static(struct lsoda_context_t * ctx, int meth);
void     scaleh(struct lsoda_context_t * ctx, int neq, double *rh, double *pdh, double hmxi);

#ifdef CFODE_STATIC
	#define cfode cfode_static
#endif

