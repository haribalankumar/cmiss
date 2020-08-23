/*
 * SuperLU interface routines
 */

#include <stdlib.h>
#include <stdio.h>

#include "dsp_defs.h"
#include "util.h"

#define SuperLU_init           superlu_init_
#define SuperLU_resetp         superlu_resetp_
#define SuperLU_destroy        superlu_destroy_
#define SuperLU_dgstrf         superlu_dgstrf_
#define SuperLU_dgstrf_refact  superlu_dgstrf_refact_
#define SuperLU_dgstrs         superlu_dgstrs_
#define Xerbla                 xerbla_

/* Local data types */
typedef struct {
  int          *etree;
  int          *perm_r;
  int          *perm_c;
  SuperMatrix  *AC;
  SuperMatrix  *L;
  SuperMatrix  *U;
} SLUSymbolic;


/* Static data */
static int local_param[10];


/* Local Routines */

SLUSymbolic *Create_SLUSymbolic( int n )
{
  SLUSymbolic *Symbolic;

  Symbolic = (SLUSymbolic *) SUPERLU_MALLOC( sizeof(SLUSymbolic) );
  if ( Symbolic == NULL ) {
    fprintf(stderr, "SUPERLU_MALLOC fails for Symbolic\n");
    exit(-1);
  }

  Symbolic->AC = SUPERLU_MALLOC( sizeof(SuperMatrix) );
  if ( Symbolic->AC == NULL ) {
    fprintf(stderr,"Malloc fails for Symbolic->AC\n");
    exit(-1);
  }
  Symbolic->L = SUPERLU_MALLOC( sizeof(SuperMatrix) );
  if ( Symbolic->L == NULL ) {
    fprintf(stderr,"Malloc fails for Symbolic->L\n");
    exit(-1);
  }
  Symbolic->U = SUPERLU_MALLOC( sizeof(SuperMatrix) );
  if ( Symbolic->U == NULL ) {
    fprintf(stderr,"Malloc fails for Symbolic->U\n");
    exit(-1);
  }

  Symbolic->etree = intMalloc(n);
  if ( Symbolic->etree == NULL ) {
    fprintf(stderr,"Malloc fails for Symbolic->etree\n");
    exit(-1);
  }
  Symbolic->perm_r = intMalloc(n);
  if ( Symbolic->perm_r == NULL ) {
    fprintf(stderr,"Malloc fails for Symbolic->perm_r\n");
    exit(-1);
  }
  Symbolic->perm_c = intMalloc(n);
  if ( Symbolic->perm_c == NULL ) {
    fprintf(stderr,"Malloc fails for Symbolic->perm_c\n");
    exit(-1);
  }

  return Symbolic;
}

void Destroy_SLUSymbolic( SLUSymbolic *Symbolic )
{
  Destroy_CompCol_Permuted( Symbolic->AC );
  Destroy_SuperNode_Matrix( Symbolic->L );
  Destroy_CompCol_Matrix( Symbolic->U );
  SUPERLU_FREE( Symbolic->AC );
  SUPERLU_FREE( Symbolic->L );
  SUPERLU_FREE( Symbolic->U );

  SUPERLU_FREE( Symbolic->etree );
  SUPERLU_FREE( Symbolic->perm_r );
  SUPERLU_FREE( Symbolic->perm_c );

  SUPERLU_FREE( Symbolic );

  return;
}


/* Public interface */

int sp_ienv(int ispec)
{
  int ierr;

  if (ispec < 1 || ispec > 6)
  {
    ierr = 1;
    Xerbla("sp_ienv", &ierr);
  }

  return local_param[ispec-1];
}


void SuperLU_resetp(int param[10])
{
  param[0] =       10; /* Panel size */
  param[1] =        5; /* Relaxation factor */
  param[2] =      100; /* Max supernode size */
  param[3] =      200; /* Min row dim */
  param[4] =       40; /* Min col dim */
  param[5] =       20; /* Fill of L, U */

  param[6] = param[5]; /* Fill of U (only used in superLU_mt)*/
  param[7] = param[5]; /* Fill of L supernodes (only used in superLU_mt)*/

  param[8] =        1; /* Column ordering */
  param[9] =        1; /* Ncpu (only used in superLU_mt)*/

  return;
}


void SuperLU_init(int param[10])
{
  int i;

  for (i = 0; i < 10; i++)
  {
    local_param[i] = abs(param[i]);
  }

  return;
}


void SuperLU_destroy(
  SLUSymbolic   **Symbolic )
{
  Destroy_SLUSymbolic( *Symbolic );
  *Symbolic = NULL;

  return;
}


void SuperLU_dgstrf(
  int           *n,
  int           *nza,
  double        *avalues,
  int           *colptr,
  int           *rowind,
  SLUSymbolic   **Symbolic,
  int           param[10],
  int           *info )
{
  int           i, lwork = 0, column_ordering, panel_size, relax;
  double        diag_pivot_thresh = 1.0, drop_tol = 0.0;
  char          *refact;
  SuperMatrix   A;
  int           *perm_c, *perm_r, *etree;
  SuperMatrix   *AC, *L, *U;

  /* Set the parameters */
  SuperLU_init(param);
  panel_size      = param[0];
  relax           = param[1];
  column_ordering = param[8];
  refact          = "N";

  /* Adjust to 0-based indexing */
  for (i = 0; i < *nza; ++i) --rowind[i];
  for (i = 0; i <=  *n; ++i) --colptr[i];

	/* Start the statistics collection */
  StatInit(panel_size, relax);

  /* Create A Super Matrix */
  dCreate_CompCol_Matrix(&A, *n, *n, *nza, avalues, rowind, colptr, NC, _D, GE);

	/* Create the sparsity information structure */
  *Symbolic = Create_SLUSymbolic( *n );

  /* Alias pointers */
  AC     = (*Symbolic)->AC;
  etree  = (*Symbolic)->etree;
  perm_c = (*Symbolic)->perm_c;
  perm_r = (*Symbolic)->perm_r;
  L      = (*Symbolic)->L;
  U      = (*Symbolic)->U;

  /* Get column permutation vector */
  get_perm_c(column_ordering, &A, perm_c);

  /* Preorder the matrix and calculate the elimation tree */
  sp_preorder(refact, &A, perm_c, etree, AC);

  /* Factorise the matrix */
  dgstrf(refact, AC, diag_pivot_thresh, drop_tol, relax, panel_size, etree,
    NULL, lwork, perm_r, perm_c, L, U, info);

  /* Clean up */
  Destroy_SuperMatrix_Store(&A);
	StatFree();

  /* Restore to 1-based indexing */
  for (i = 0; i < *nza; ++i) ++rowind[i];
  for (i = 0; i <=  *n; ++i) ++colptr[i];

  return;
}


void SuperLU_dgstrf_refact(
  int           *n,
  int           *nza,
  double        *avalues,
  int           *colptr,
  int           *rowind,
  SLUSymbolic   **Symbolic,
  int           param[10],
  int           *info )
{
  int           i, lwork = 0, panel_size, relax;
  double        diag_pivot_thresh = 1.0, drop_tol = 0.0;
  int           *perm_c, *perm_r, *etree;
  SuperMatrix   *AC, *L, *U;
  NCformat      *ACStore;
  char          *refact;

  /* Set the parameters */
  panel_size      = param[0];
  relax           = param[1];
  refact          = "Y";

  /* Alias pointers */
  AC     = (*Symbolic)->AC;
  etree  = (*Symbolic)->etree;
  perm_c = (*Symbolic)->perm_c;
  perm_r = (*Symbolic)->perm_r;
  L      = (*Symbolic)->L;
  U      = (*Symbolic)->U;

  /* Update the permutation values */
  ACStore         = AC->Store;
  ACStore->nzval  = avalues;

  /* Adjust to 0-based indexing */
  for (i = 0; i < *nza; ++i) --rowind[i];
  for (i = 0; i <=  *n; ++i) --colptr[i];

	/* Start the statistics collection */
  StatInit(panel_size, relax);

  /* Factorise the matrix */
  dgstrf(refact, AC, diag_pivot_thresh, drop_tol, relax, panel_size, etree,
    (void *) NULL, lwork, perm_r, perm_c, L, U, info);

	/* Clean Up */
	StatFree();

  /* Restore to 1-based indexing */
  for (i = 0; i < *nza; ++i) ++rowind[i];
  for (i = 0; i <=  *n; ++i) ++colptr[i];

  return;
}


void SuperLU_dgstrs(
  char          *trans,
  int           *n,
  double        *bvalues,
  int           *nrhs,
  int           *ldb,
  SLUSymbolic   **Symbolic,
  int           param[10],
  int           *info,
  long          trans_len )
{
  SuperMatrix   B;
	int           panel_size, relax;
  char          *transp;
  int           *perm_c, *perm_r;
  SuperMatrix   *L, *U;

  /* Set the parameters */
  panel_size = param[0];
  relax      = param[1];

  /* Alias pointers */
  perm_c = (*Symbolic)->perm_c;
  perm_r = (*Symbolic)->perm_r;
  L      = (*Symbolic)->L;
  U      = (*Symbolic)->U;

  /* Check value of trans */
  if (trans_len < 1) {
    transp = "t";
  }
  else {
    transp = trans;
  }

	/* Start the statistics collection */
  StatInit(panel_size, relax);

  /* Create B Super Matrix */
  dCreate_Dense_Matrix(&B, *n, *nrhs, bvalues, *ldb, DN, _D, GE);

  /* Solve */
  dgstrs(transp, L, U, perm_r, perm_c, &B, info);

  /* Clean up and go */
  Destroy_SuperMatrix_Store(&B);
	StatFree();

  return;
}
