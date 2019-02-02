
#ifndef __ITSOL_ILUK_H__
#define __ITSOL_ILUK_H__

#include "type-defs.h"

int lofC( int lofM, csptr csmat, iluptr lu, FILE *fp ); 
int ilukC( int lofM, csptr csmat, iluptr lu, FILE *fp );

#endif
