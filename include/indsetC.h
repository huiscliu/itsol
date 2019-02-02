
#ifndef __ITSOL_INDSETC_H__
#define __ITSOL_INDSETC_H__

#include "sets.h"
#include "misc.h"

int add2is(int *last, int nod, int *iord, int *riord);
int add2com(int *nback, int nod, int *iord, int *riord);
int indsetC(csptr mat, int bsize, int *iord, int *nnod, double tol) ;
int weightsC(csptr mat, double *w);

#endif
