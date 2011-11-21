#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lsoda.h"

void lsoda_(void (*)(int *, double *, double *, double *),
		     int *, double *, double *, double *,
		     int *, double *, double *, int *, int *,
		     int *, double *,int *,int *, int *,
		     void (*)(int *, double *, double *, int *,
			      int *, double *, int *),
		     int *);

void rwarnc_(char * str, int *len) {
	int l = *len;
	while(--l) putchar(*(str++));
}

/**
 * This is a c wrapper of lsoda.f
 * tcrit is the critical time that the integrator shall stop, can be NULL
 * jac is the jacobian function, can be NULL
 * hmin is the min step can be NULL
 * hmax is the max step can be NULL
 *
 * y is the init condition, yout is the final value
 * t is the initial time, overwritten with the ending time
 * tout is the final time.
 * derives is the function calculating ydot
 * rtol[lrtol] is the array for relative error tolarance. can be length 1 if the rtol is same for all dimensions.
 * atol[larol] .................aboslute ....................................
 * 
 * returns the state code of losda
 *
 * this function is adapted from source code of R package odesolve.
 * 
 * This code is in public domain
 * Yu Feng 2011, Carnegie Mellon Unviersity
 */

int lsoda(double * y, double * yout, int neq, double * t, double tout, deriv_func *derivs, double * rtol, int lrtol,
		double * atol, int latol, double * tcrit, jac_func * jac, double * hmin, double * hmax)
{
  int i, j, k, repcount, istate;
  double  *rwork, tin;
  int itol, itask, iopt, lrw, liw, *iwork, jt, lrn, lrs;

  for (j = 0; j < neq; j++) yout[j] = y[j];

  if (latol == 1 && lrtol == 1 ) itol = 1;
  if (latol  > 1 && lrtol == 1 ) itol = 2;
  if (latol == 1 && lrtol  > 1 ) itol = 3;
  if (latol  > 1 && lrtol  > 1 ) itol = 4;
  if (tcrit != NULL)
    {
      itask = 4;
    }
  else
    {
      itask = 1;
    }
  istate = 1;
  iopt = 0;
  lrn = 20 + 16 * neq;
  lrs = 22 + 9 * neq + neq * neq;
  if (lrn > lrs) lrw = lrn;
  else lrw = lrs;

  rwork = (double*) alloca(lrw * sizeof(double));
  memset(rwork, 0, sizeof(double) * lrw);

  if (itask == 4) rwork[0] = *tcrit;
  rwork[5] = hmax?*hmax:0;
  rwork[6] = hmin?*hmin:0;
  if (rwork[5] >0 || rwork[6] >0) iopt = 1;

  liw = 20 + neq;

  iwork = (int *) alloca(liw * sizeof(int));
  memset(iwork, 0, sizeof(int) * liw);

  if (jac!= NULL)
    jt = 1;
  else
    jt = 2;

      /* I still need to trap possible error returns of lsoda */
      /* based on return values in istate, */
      tin = *t;
      repcount = 0;

      do
	{
	  if (istate == -1) istate = 2;
	  if (istate == -2)
	    {
	      for (j = 0; j < lrtol; j++) rtol[j] *= 10.0;
	      for (j = 0; j < latol; j++) atol[j] *= 10.0;
	      fprintf(stderr, "Excessive precision requested.  `rtol' and `atol' have been scaled upwards by the factor %g\n",10.0);
	      istate = 3;
	      
	    }
	  lsoda_ (derivs, &neq, yout, &tin, &tout,
			   &itol, rtol, atol, &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jac, &jt);
	  repcount ++;
	} while (tin < tout && istate >= -2 && repcount < 20);

      if (istate == -3)
	{
	  error("Illegal input to lsoda\n");
	}
      else
	{
	  *t = tin;
	}

  return istate;
}


