int stoda(struct common_t * common, int neq, double *y, _lsoda_f f, void *_data, int jstart, double hmxi, double hmin, int mxords, int msordn);
int      correction(struct common_t * common, int neq, double *y, _lsoda_f f, double pnorm, double *del, double *delp, double *told,
						   int *ncf, double *rh, int *m, double hmin, void *_data);
int prja(struct common_t * common, int neq, double *y, _lsoda_f f, void *_data);

int      ewset(struct common_t * common, const int neq, double ewt[], const double rtol[], const double atol[], const double *ycur);
void     solsy(struct common_t * common, int neq, double *y);
int      orderswitch(struct common_t * common, int neq, double rhup, double dsm, double *pdh, double *rh, int kflag, int maxord);
void     intdy(struct common_t * common, int neq, double t, int k, double *dky, int *iflag);
int      corfailure(struct common_t * common, int neq, double *told, double *rh, int *ncf, double hmin);
void     methodswitch(struct common_t * common, int neq, double dsm, double pnorm, double *pdh, double *rh, int mxords, int mxordn);
void     cfode(struct common_t * common, int meth);
void     cfode_static(struct common_t * common, int meth);
void     scaleh(struct common_t * common, int neq, double *rh, double *pdh, double hmxi);

#ifdef CFODE_STATIC
	#define cfode cfode_static
#endif

