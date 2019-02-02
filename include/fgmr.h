
#ifndef __ITSOL_FGMR_H__
#define __ITSOL_FGMR_H__

#include "type-defs.h"

/*----------------------------------------------------------------------
|                 *** Preconditioned FGMRES ***                  
+-----------------------------------------------------------------------
| This is a simple version of the ARMS preconditioned FGMRES algorithm. 
+-----------------------------------------------------------------------
| Y. S. Dec. 2000. -- Apr. 2008   -- Jul. 2009 
+-----------------------------------------------------------------------
| on entry:
|---------- 
|
|(Amat)   = matrix struct. the matvec operation is Amat->matvec.
|(lu)     = preconditioner struct.. the preconditioner is lu->precon
|           if (lu == NULL) the no-preconditioning option is invoked.
| rhs     = real vector of length n containing the right hand side.
| sol     = real vector of length n containing an initial guess to the
|           solution on input.
| tol     = tolerance for stopping iteration
| im      = Krylov subspace dimension 
| (itmax) = max number of iterations allowed. 
| fits    = NULL: no output
|        != NULL: file handle to output " resid vs time and its" 
|
| on return:
|---------- 
| fgmr      int =  0 --> successful return.
|           int =  1 --> convergence not achieved in itmax iterations.
| sol     = contains an approximate solution (upon successful return).
| itmax   = has changed. It now contains the number of steps required
|           to converge -- 
+-----------------------------------------------------------------------
| internal work arrays:
|----------       
| vv      = work array of length [im+1][n] (used to store the Arnoldi
|           basis)
| hh      = work array of length [im][im+1] (Arnoldi matrix)
| z       = work array of length [im][n] to store preconditioned vectors
+-----------------------------------------------------------------------
| subroutines called :
|     matvec and
|     preconditionning operation 
+---------------------------------------------------------------------*/
int fgmr(SMatptr Amat, SPreptr lu, double *rhs, double *sol, 
         double tol, int im, int *itmax, FILE *fits);

#endif
