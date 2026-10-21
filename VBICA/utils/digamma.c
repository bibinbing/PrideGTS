/*
 * digamma.c
 *
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"

extern double psi_(double *xx);

#define Xin prhs[0]

void Psi (double y[], double x[])
{
  y[0]=psi_(&x[0]);
  /*y[0]=2.0*x[0];*/
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*
   * Gateway function to matlab.
   *
   * y = digamma(x)
   */

  double *x,*y;
  int m,n;
  
  if (nrhs < 1)
    mexErrMsgTxt("digamma requires at least one input argument.");
  if (nlhs < 1)
    mexErrMsgTxt("digamma requires at least one output argument.");

  m = mxGetM(Xin);
  n = mxGetN(Xin);
  if (!mxIsNumeric(Xin) || mxIsComplex(Xin) || 
      !mxIsDouble(Xin) || !(m== 1 && n==1)) 
    mexErrMsgTxt("digamma needs a noncomplex scale double as input.");

  /* Create matrix for the return argument */
  plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);

  /* Assign pointers to input and output */
  x = mxGetPr(Xin);
  y = mxGetPr(plhs[0]);
  
  Psi(y,x);
  /*printf("x=%1.2f, y=%1.2f",*x,*y);*/
  return;

}

  
