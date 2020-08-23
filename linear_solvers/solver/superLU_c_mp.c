/*
 * SuperLU_mt interface routines
 */

#include <stdlib.h>
#include <stdio.h>

#include "pdsp_defs.h"
#include "util.h"
#include "machines.h"

#define SuperLU_init           superlu_init_
#define SuperLU_resetp         superlu_resetp_
#define SuperLU_destroy        superlu_destroy_
#define SuperLU_dgstrf         superlu_dgstrf_
#define SuperLU_dgstrf_refact  superlu_dgstrf_refact_
#define SuperLU_dgstrs         superlu_dgstrs_
#define Xerbla                 xerbla_

/* Local data types */

typedef struct {
  int               *perm_r;
  int               *perm_c;
	SuperMatrix       *AC;
	SuperMatrix       *L;
	SuperMatrix       *U;
  pdgstrf_options_t *pdgstrf_options;
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

  Symbolic->pdgstrf_options = (pdgstrf_options_t *)
    SUPERLU_MALLOC( sizeof(pdgstrf_options_t) );
  if ( Symbolic->pdgstrf_options == NULL ) {
    fprintf(stderr, "SUPERLU_MALLOC fails for Symbolic->pdgstrf_options\n");
    exit(-1);
  }

  return Symbolic;
}

void Destroy_SLUSymbolic( SLUSymbolic *Symbolic )
{
  Destroy_CompCol_Permuted( Symbolic->AC );
  Destroy_SuperNode_SCP( Symbolic->L );
  Destroy_CompCol_NCP( Symbolic->U );
  SUPERLU_FREE( Symbolic->AC );
  SUPERLU_FREE( Symbolic->L );
  SUPERLU_FREE( Symbolic->U );

  SUPERLU_FREE( Symbolic->perm_r );
  SUPERLU_FREE( Symbolic->perm_c );

  SUPERLU_FREE( Symbolic->pdgstrf_options->etree );
  SUPERLU_FREE( Symbolic->pdgstrf_options->colcnt_h );
  SUPERLU_FREE( Symbolic->pdgstrf_options->part_super_h );
  SUPERLU_FREE( Symbolic->pdgstrf_options );

  SUPERLU_FREE( Symbolic );

  return;
}

/* Public interface */

int sp_ienv(int ispec)
{
  int ierr;

  if (ispec < 1 || ispec > 8)
  {
    ierr = 1;
    Xerbla("sp_ienv", &ierr);
  }

  return local_param[ispec-1];
}


void SuperLU_resetp(int param[10])
{
#if ( MACH==SGI )
  param[0] =  20; /* Panel size */
  param[1] =   1; /* Relaxation factor */
  param[2] = 100; /* Max supernode size */
  param[3] = 800; /* Min row dim */
  param[4] = 100; /* Min col dim */
#elif ( MACH==ORIGIN ) 
  param[0] =  12;
  param[1] =   1;
  param[2] = 100;
  param[3] = 400;
  param[4] = 100;
#elif ( MACH==DEC )
  param[0] =  16;
  param[1] =   1;
  param[2] =  50;
  param[3] = 100;
  param[4] =  40;
#elif ( MACH==CRAY_PVP )
  param[0] =   1;
  param[1] =   1;
  param[2] =  64;
  param[3] = 400;
  param[4] = 200;
#elif ( MACH==SUN )
  param[0] =   8;
  param[1] =   1;
  param[2] = 100;
  param[3] = 400;
  param[4] =  40;
#elif ( MACH==IBM )
  param[0] =  20;
  param[1] =   1;
  param[2] = 100;
  param[3] = 800;
  param[4] = 100;
#else
  param[0] =   8; /* Panel size */
  param[1] =   1; /* Relaxation factor */
  param[2] = 100; /* Max supernode size */
  param[3] = 200; /* Min row dim */
  param[4] =  40; /* Min col dim */
#endif
  param[5] = -20; /* Fill of L */
  param[6] = -20; /* Fill of U */
  param[7] = -10; /* Fill of L supernodes */

  param[8] =   1; /* Column ordering */
  param[9] =   1; /* Ncpu */

  return;
}


void SuperLU_init(int param[10])
{
  int i;

  for (i = 0; i <  5; i++) {
    local_param[i] =  param[i];
  }
	for (i = 5; i <  8; i++) {
		local_param[i] = -param[i];
	}
	for (i = 8; i < 10; i++) {
		local_param[i] =  param[i];
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
  int               i, lwork = 0;
  int               column_ordering, panel_size, relax, num_proc;
  double            diag_pivot_thresh = 1.0, drop_tol = 0.0;
  yes_no_t          refact = NO, usepr=NO;
  pdgstrf_options_t *pdgstrf_options;
  SuperMatrix       A;
  int               *perm_c, *perm_r;
  SuperMatrix       *AC, *L, *U;
	Gstat_t           Gstat;


  /* Set the parameters */
  SuperLU_init(param);
  panel_size      = param[0];
  relax           = param[1];
  column_ordering = param[8];
  num_proc        = param[9];
  refact          = NO;

  /* Adjust to 0-based indexing */
  for (i = 0; i < *nza; ++i) --rowind[i];
  for (i = 0; i <=  *n; ++i) --colptr[i];

  /* Create A Super Matrix */
  dCreate_CompCol_Matrix(&A, *n, *n, *nza, avalues, rowind, colptr, NC, _D, GE);

	/* Create the sparsity information structure */
  *Symbolic = Create_SLUSymbolic( *n );

  /* Alias pointers */
  AC     = (*Symbolic)->AC;
  perm_c = (*Symbolic)->perm_c;
  perm_r = (*Symbolic)->perm_r;
  L      = (*Symbolic)->L;
  U      = (*Symbolic)->U;
	pdgstrf_options = (*Symbolic)->pdgstrf_options;

	/* Start the statistics collection */
  StatAlloc(*n, num_proc, panel_size, relax, &Gstat);
  StatInit(*n, num_proc, &Gstat);

  /* Get column permutation vector */
  get_perm_c(column_ordering, &A, perm_c);

  /* Initialise the option structure pdgstrf_options using the
     user-input parameters */
  pdgstrf_init(num_proc, refact, panel_size, relax,
               diag_pivot_thresh, usepr, drop_tol, perm_c, perm_r,
               NULL, lwork, &A, AC, pdgstrf_options, &Gstat);

  /* Factorise the matrix */
  pdgstrf(pdgstrf_options, AC, perm_r, L, U, &Gstat, info);

  /* Clean up */
  Destroy_SuperMatrix_Store(&A);
	StatFree(&Gstat);

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
  int           i;
  int           panel_size, relax, num_proc;
  pdgstrf_options_t *pdgstrf_options;
  int           *perm_r;
  SuperMatrix   *AC, *L, *U;
	Gstat_t       Gstat;


  /* Set the parameters */
  panel_size      = param[0];
  relax           = param[1];
  num_proc        = param[9];

  /* Adjust to 0-based indexing */
  for (i = 0; i < *nza; ++i) --rowind[i];
  for (i = 0; i <=  *n; ++i) --colptr[i];

  /* Alias pointers */
  AC     = (*Symbolic)->AC;
  perm_r = (*Symbolic)->perm_r;
  L      = (*Symbolic)->L;
  U      = (*Symbolic)->U;
	pdgstrf_options = (*Symbolic)->pdgstrf_options;
	pdgstrf_options->refact = YES;

	/* Start the statistics collection */
  StatAlloc(*n, num_proc, panel_size, relax, &Gstat);
  StatInit(*n, num_proc, &Gstat);

  /* Factorise the matrix */
  pdgstrf(pdgstrf_options, AC, perm_r, L, U, &Gstat, info);

	/* Clean Up */
	StatFree(&Gstat);

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
  int           panel_size, relax, num_proc;
  trans_t       transpose;
  int           *perm_c, *perm_r;
  SuperMatrix   *L, *U;
  Gstat_t       Gstat;

  /* Set the parameters */
  panel_size = param[0];
  relax      = param[1];
  num_proc   = param[9];

  /* Alias pointers */
  perm_c = (*Symbolic)->perm_c;
  perm_r = (*Symbolic)->perm_r;
  L      = (*Symbolic)->L;
  U      = (*Symbolic)->U;

  /* Set the transpose flag */
  if (trans_len < 1) {
    transpose = TRANS;
  }
  else if (trans[0] == 'N' || trans[0] == 'n') {
    transpose = NOTRANS;
  }
  else {
    transpose = TRANS;
  }

  /* Start the statistics collection */
  StatAlloc(*n, num_proc, panel_size, relax, &Gstat);
  StatInit(*n, num_proc, &Gstat);

  /* Create B Super Matrix */
  dCreate_Dense_Matrix(&B, *n, *nrhs, bvalues, *ldb, DN, _D, GE);

  /* Solve */
  dgstrs(transpose, L, U, perm_r, perm_c, &B, &Gstat, info);

  /* Clean up and go */
  Destroy_SuperMatrix_Store(&B);
  StatFree(&Gstat);

  return;
}
