struct vec_t {
	double **yh, **wm, *ewt, *savf, *acor;
	int     *ipvt;
};

int stoda(int neq, double *y, _lsoda_f f, void *_data, int jstart, double hmxi, double hmin, int mxords, int msordn);
void     correction(int neq, double *y, _lsoda_f f, int *corflag, double pnorm, double *del, double *delp, double *told,
						   int *ncf, double *rh, int *m, double hmin, void *_data);
void     prja(int neq, double *y, _lsoda_f f, void *_data);

int ewset(double ewt[], int itol, double rtol[], double atol[], double *ycur, int n);
void     resetcoeff(void);
void     solsy(double *y, double ** wm, int * ipvt, int n);
int      orderswitch(double *rhup, double dsm, double *pdh, double *rh, int kflag);
void     intdy(double t, int k, double *dky, int *iflag);
void     corfailure(double *told, double *rh, int *ncf, int *corflag, double hmin);
void     methodswitch(double dsm, double pnorm, double *pdh, double *rh, int mxords, int mxordn);
void     cfode(int meth);
void     scaleh(double *rh, double *pdh, double hmxi);


