
#ifndef __ITSOL_AUXILL_H__
#define __ITSOL_AUXILL_H__

#include "type-defs.h"
#include "sets.h"

int read_inputs(char *in_file, io_t *pio);

int get_matrix_info(FILE *fmat, io_t *pio);
void output_blocks(int nBlock, int *nB, FILE *f);
void output_perm(int n, int *perm, FILE *f);
int read_coo(double **VAL, int **COL, int **ROW, io_t *pio, double **rhs, double **sol, int job);

int readhb_c(int *NN, double **AA, int **JA, int **IA, io_t *pio, 
	     double **rhs, double **sol, int *rsa);

/*-----------------------------------------------------------*
   Output AA, JA, IA
   fmt == 0, output in CSC
   fmt == 1, output in CSR
 *-----------------------------------------------------------*/
int readhb_2(int *NN, double **AA, int **JA, int **IA, io_t *pio, double **rhs, double **sol, int *rsa, int fmt);

void output_header(io_t *pio);

void output_header_vb(io_t *pio);

void output_result(int lfil, io_t *pio, int iparam);

/*-------------------------------------------------*/
/* sets parameters required by arms preconditioner */
/* input io_t, Dscale                            */
/* output ipar tolcoef, lfil                       */
/*-------------------- trigger an error if not set */
void set_arms_pars(io_t* io, int Dscale, int *ipar, double *dropcoef, int *lfil);

void randvec (double *v, int n);

#endif
