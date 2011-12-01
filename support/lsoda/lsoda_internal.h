int stoda(struct lsoda_context_t * ctx, double *y, int jstart);
int      correction(struct lsoda_context_t * ctx, double *y, double pnorm, double *del, double *delp, double *told,
						   int *ncf, double *rh, int *m, double hmin);
int prja(struct lsoda_context_t * ctx, double *y);

int      ewset(struct lsoda_context_t * ctx, double ewt[], const double rtol[], const double atol[], const double *ycur);
void     solsy(struct lsoda_context_t * ctx, double *y);
int      orderswitch(struct lsoda_context_t * ctx, double rhup, double dsm, double *pdh, double *rh, int kflag, int maxord);
int      intdy(struct lsoda_context_t * ctx, double t, int k, double *dky);
int      corfailure(struct lsoda_context_t * ctx, double *told, double *rh, int *ncf, double hmin);
void     methodswitch(struct lsoda_context_t * ctx, double dsm, double pnorm, double *pdh, double *rh);
void     cfode(struct lsoda_context_t * ctx, int meth);
void     cfode_static(struct lsoda_context_t * ctx, int meth);
void     scaleh(struct lsoda_context_t * ctx, double *rh, double *pdh);

#ifdef CFODE_STATIC
	#define cfode cfode_static
#endif

