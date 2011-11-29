struct vec_t {
	double **yh, **wm, *ewt, *savf, *acor;
	int     *ipvt;
};

int stoda(int neq, double *y, _lsoda_f f, void *_data, int jstart, double hmxi, double hmin, int mxords, int msordn);
int      correction(int neq, double *y, _lsoda_f f, double pnorm, double *del, double *delp, double *told,
						   int *ncf, double *rh, int *m, double hmin, void *_data);
int prja(int neq, double *y, _lsoda_f f, void *_data);

int      ewset(int neq, double ewt[], int itol, double rtol[], double atol[], double *ycur);
void     resetcoeff(void);
void     solsy(int neq, double *y, double ** wm, int * ipvt);
int      orderswitch(int neq, double rhup, double dsm, double *pdh, double *rh, int kflag);
void     intdy(int neq, double t, int k, double *dky, int *iflag);
int      corfailure(int neq, double *told, double *rh, int *ncf, double hmin);
void     methodswitch(int neq, double dsm, double pnorm, double *pdh, double *rh, int mxords, int mxordn);
void     cfode(int meth);
void     scaleh(int neq, double *rh, double *pdh, double hmxi);


