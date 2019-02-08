
#include "auxill.h"

#define ERR_AUXIL  10
#define ITS_MAX_NUM_LEV 10          /* maximum number of levels for arms    */

/*-------------------- reads a matrix in coordinate format. 
  !  arrays VAL, COL, ROW are allocated and created 
  !  for rhs: memory allocation done + artificial rhs created.
  !  various other things are filled in pio  
  ! job = 0  - want C indexing 
  ! job = 1  - want FORTRAN indexing 
  !------------------------------------------------------------*/
int itsol_read_coo(double **VAL, int **COL, int **ROW, double **rhs, double **sol, char *Fname, int *ndim, int *nnnz)
{
    FILE *matf = NULL;
    double *aa;
    int *ii, *jj;
    int k, n, nnz;
    char str[ITS_MAX_LINE];
    /*-------------------- start */
    if ((matf = fopen(Fname, "r")) == NULL) {
        fprintf(stdout, "Cannot Open Matrix\n");
        return (ERR_AUXIL + 3);
    }
    /*-------------------- mtx format .. in some cases n, 
      nnz are read separately per line */
    /*-------------------- try a 100 lines max of comments */
    for (k = 0; k < 100; k++) {
        fgets(str, ITS_MAX_LINE, matf);
        if (memcmp(str, "%", sizeof(char)) != 0)
            break;
    }
    if (k == 99)
        return (ERR_AUXIL + 3);
    sscanf(str, " %d %d %d", &n, &k, &nnz);
    if (n != k) {
        fprintf(stdout, "This is not a square matrix -- stopping \n");
        return (ERR_AUXIL + 4);
    }
    /* separate reads for n and nnz 
       fscanf(matf," %d", &n); 
       fscanf(matf," %d", &nnz);  
     */
    *ndim = n;
    *nnnz = nnz;

    /*-------------------- allocate memory for matrix and rhs --- */
    *rhs = (double *)itsol_malloc(n * sizeof(double), "read_coo:1");
    *sol = (double *)itsol_malloc(n * sizeof(double), "read_coo:2");
    aa = (double *)itsol_malloc(nnz * sizeof(double), "read_coo:3");
    jj = (int *)itsol_malloc(nnz * sizeof(int), "read_coo:4");
    ii = (int *)itsol_malloc(nnz * sizeof(int), "read_coo:5");
    /*-------------------- long live fortran77 --- */
    for (k = 0; k < nnz; k++) {
        fscanf(matf, "%d  %d  %s", &ii[k], &jj[k], str);
        aa[k] = atof(str);
    }

    *ROW = ii;
    *COL = jj;
    *VAL = aa;
    /*-------------------- TO UPDATE 
      if rhs and sols are available load them -- otherwise
      generate artificial ones. */
    for (k = 0; k < n; k++) {
        (*sol)[k] = 1.0;
        (*rhs)[k] = 0.0;
    }
    for (k = 0; k < nnz; k++)
        (*rhs)[ii[k]] += aa[k] * (*sol)[jj[k]];

    fclose(matf);
    return (0);
}
