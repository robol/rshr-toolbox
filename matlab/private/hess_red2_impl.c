#include "mex.h"
#include <string.h>

extern void hr_impl_full_to_kh_ (double * A, long int * n, double * U, 
				 double * V, long int * k);

extern void hr_impl_kh_to_h_ (double * A, long int * n, double * U, 
			      double * V, long int * k, double * dd, 
			      double * ss);

void mexFunction (int nlhs, mxArray * plhs[], 
		  int nrhs, const mxArray * prhs[])
{
  /* No need to accurately check dimensions, they have already
   * been checked in the parent function. */
  long int n = mxGetM (prhs[1]);
  long int k = mxGetN (prhs[1]);

  plhs[0] = mxCreateDoubleMatrix (n, 2 * k + 3, mxREAL);
  plhs[1] = mxCreateDoubleMatrix (n, k, mxREAL);
  plhs[2] = mxCreateDoubleMatrix (n, k, mxREAL);
  plhs[3] = mxCreateDoubleMatrix (n, 1, mxREAL);
  plhs[4] = mxCreateDoubleMatrix (n-1,1, mxREAL);

  double * A = mxGetPr (plhs[0]);
  double * U = mxGetPr (plhs[1]);
  double * V = mxGetPr (plhs[2]);
  double * dd = mxGetPr (plhs[3]);
  double * ss = mxGetPr (plhs[4]);

  memcpy (A, mxGetPr (prhs[0]), n * (2 * k + 3) * sizeof (double));
  memcpy (U, mxGetPr (prhs[1]), n * k * sizeof (double));
  memcpy (V, mxGetPr (prhs[2]), n * k * sizeof (double)); 

  /* Perform reduction from diagonal plus low rank to k-Hessenberg
   * form as k-banded plus low rank with U = [ X ; 0 ]. */
  hr_impl_full_to_kh_ (A, &n, U, V, &k);

  /* Complete Hessenberg reduction by taking the matrix to standard 
   * Hessenberg form */
  hr_impl_kh_to_h_ (A, &n, U, V, &k, dd, ss);
}
