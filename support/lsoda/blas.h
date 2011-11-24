
double ddot(int n, double dx[], int incx, double dy[], int incy);
void dgesl(double **a, int n, int * ipvt, double b[], int job);
void dgefa(double ** a, int n, int * ipvt, int * info);
void daxpy(int n, double da, double dx[], int incx, double dy[], int incy);
int idamax(int n, double dx[], int incx);
void dscal(int n, double da, double dx[], int incx);
double vmnorm(int n, double *v, double *w);
double fnorm(int n, double **a, double *w);
