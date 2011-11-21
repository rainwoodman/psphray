typedef void deriv_func(int * neq, double * t, double * y,double * ydot);
typedef void jac_func(int * neq, double * t, double * y, int * ml,
		      int * mu, double * pd, int * nrowpd);

int lsoda(double * y, double * yout, int neq, double * t, double tout, deriv_func *derivs, double * rtol, int lrtol,
		double * atol, int latol, double * tcrit, jac_func * jac, double * hmin, double * hmax);

