
#ifndef __ITSOL_AUXILL_H__
#define __ITSOL_AUXILL_H__

#include "type-defs.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

int itsol_read_inputs(char *in_file, ITS_IOT *pio);

int itsol_get_matrix_info(FILE *fmat, ITS_IOT *pio);
void itsol_output_blocks(int nBlock, int *nB, FILE *f);
void itsol_output_perm(int n, int *perm, FILE *f);
int itsol_read_coo(double **VAL, int **COL, int **ROW, ITS_IOT *pio, double **rhs, double **sol, int job);

int itsol_readhb_c(int *NN, double **AA, int **JA, int **IA, ITS_IOT *pio, 
	     double **rhs, double **sol, int *rsa);

/*-----------------------------------------------------------*
   Output AA, JA, IA
   fmt == 0, output in CSC
   fmt == 1, output in CSR
 *-----------------------------------------------------------*/
int itsol_readhb_2(int *NN, double **AA, int **JA, int **IA, ITS_IOT *pio, double **rhs, double **sol, int *rsa, int fmt);

void itsol_output_header(ITS_IOT *pio);

void itsol_output_header_vb(ITS_IOT *pio);

void itsol_output_result(int lfil, ITS_IOT *pio, int iparam);

/*-------------------------------------------------*/
/* sets parameters required by arms preconditioner */
/* input ITS_IOT, Dscale                            */
/* output ipar tolcoef, lfil                       */
/*-------------------- trigger an error if not set */
void itsol_set_arms_pars(ITS_PARS *io, int Dscale, int *ipar, double *dropcoef, int *lfil);

ITS_PARS itsol_set_iot_to_pars(ITS_PARS *io);

void itsol_randvec (double *v, int n);

#ifdef __cplusplus
}
#endif
#endif
