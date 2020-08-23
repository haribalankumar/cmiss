/*
 * UMFPACK 4.0 interface routines
 */

#include <stdlib.h>
#include <stdio.h>

#include "umfpack.h"

#define Umfpack4_Factor   umfpack4_factor_c_
#define Umfpack4_Refactor umfpack4_refactor_c_
#define Umfpack4_Solve    umfpack4_solve_c_
#define Umfpack4_Free     umfpack4_free_c_
#define Umfpack4_Error    umfpack4_error_


/* Private functions */

static void copy_out(char *to, char *from, int tolen)
{
  int i = 0;

  for (i=0;i<tolen;i++)
  {
    if (from[i] == '\0') {
      break;
    }
    to[i] = from[i];
  }

  for (i=i;i<tolen;i++)
  {
    to[i] = ' ';
  }
}

/* Public functions */

void Umfpack4_Error(int *ierr, char *error, int errlen)
{
  char  buff[32];

  switch (*ierr) {
  case UMFPACK_WARNING_singular_matrix:
    {
      copy_out(error, "Singular Matrix", errlen);
      break;
    }
  case UMFPACK_ERROR_out_of_memory:
    {
      copy_out(error, "Out of memory", errlen);
      break;
    }
  case UMFPACK_ERROR_invalid_Numeric_object:
    {
      copy_out(error, "Invalid Numeric object", errlen);
      break;
    }
  case UMFPACK_ERROR_invalid_Symbolic_object:
    {
      copy_out(error, "Invalid Symbolic object", errlen);
      break;
    }
  case UMFPACK_ERROR_argument_missing:
    {
      copy_out(error, "Argument missing", errlen);
      break;
    }
  case UMFPACK_ERROR_n_nonpositive:
    {
      copy_out(error, "N nonpositive", errlen);
      break;
    }
  case UMFPACK_ERROR_nz_negative:
    {
      copy_out(error, "NZ negative", errlen);
      break;
    }
  case UMFPACK_ERROR_jumbled_matrix:
    {
      copy_out(error, "Jumbled matrix", errlen);
      break;
    }
  case UMFPACK_ERROR_Ap0_nonzero:
    {
      copy_out(error, "Ap0 nonzero", errlen);
      break;
    }
  case UMFPACK_ERROR_row_index_out_of_bounds:
    {
      copy_out(error, "Row index out of bounds", errlen);
      break;
    }
  case UMFPACK_ERROR_different_pattern:
    {
      copy_out(error, "Different pattern", errlen);
      break;
    }
  case UMFPACK_ERROR_col_length_negative:
    {
      copy_out(error, "Col length negative", errlen);
      break;
    }
  case UMFPACK_ERROR_invalid_system:
    {
      copy_out(error, "Invalid system", errlen);
      break;
    }
  case UMFPACK_ERROR_invalid_triplet:
    {
      copy_out(error, "Invalid triplet", errlen);
      break;
    }
  case UMFPACK_ERROR_invalid_permutation:
    {
      copy_out(error, "Invalid permutation", errlen);
      break;
    }
  case UMFPACK_ERROR_problem_too_large:
    {
      copy_out(error, "Problem too large", errlen);
      break;
    }
  case UMFPACK_ERROR_internal_error:
    {
      copy_out(error, "Internal error", errlen);
      break;
    }
  default:
    {
      sprintf(buff, "Unknown error %d", *ierr);
      copy_out(error, buff, errlen);
      break;
    }
  }

  return;
}


/* Public interface */

void Umfpack4_Factor( double *A, int *n, int *nza, int *isc, int *isr,
                      void **Symbolic, void **Numeric, double *umfdef,
                      double *Control, double *Info, char *name, int *ierr,
                      int namelen )
{
  int i;

  /* Set the defaults */
  /* printf("Set the controls (for the heart of the sun/pelvis) ...\n"); */
  umfpack_di_defaults( Control );

  if (umfdef[ 0] >= 0.0) Control[UMFPACK_PIVOT_TOLERANCE]       = umfdef[ 0];
  if (umfdef[ 1] >= 0.0) Control[UMFPACK_BLOCK_SIZE]            = umfdef[ 1];
  if (umfdef[10] >= 0.0) Control[UMFPACK_DENSE_ROW]             = umfdef[10];
  if (umfdef[11] >= 0.0) Control[UMFPACK_DENSE_COL]             = umfdef[11];
  if (umfdef[12] >= 0.0) Control[UMFPACK_RELAXED_AMALGAMATION]  = umfdef[12];
  if (umfdef[13] >= 0.0) Control[UMFPACK_RELAXED2_AMALGAMATION] = umfdef[13];
  if (umfdef[13] >= 0.0) Control[UMFPACK_RELAXED3_AMALGAMATION] = umfdef[13];
  Control[UMFPACK_IRSTEP] = 0.0;

  /* Adjust to 0-based indexing */
  for (i = 0; i <=  *n; ++i) --isr[i];
  for (i = 0; i < *nza; ++i) --isc[i];

  *ierr = umfpack_di_symbolic(*n, *n, isr, isc, Symbolic, Control, Info);
  if (*ierr != UMFPACK_OK) {
    copy_out(name, "umfpack_symbolic", namelen);
    return;
  }

  *ierr = umfpack_di_numeric(isr, isc, A, *Symbolic, Numeric, Control, Info);
  if (*ierr != UMFPACK_OK) {
    copy_out(name, "umfpack_numeric", namelen);
    return;
  }

  /* Restore to 1-based indexing */
  for (i = 0; i <=  *n; ++i) ++isr[i];
  for (i = 0; i < *nza; ++i) ++isc[i];

  return;
}


void Umfpack4_Refactor( double *A, int *n, int *nza, int *isc, int *isr,
                        void **Symbolic, void **Numeric, double *Control,
                        double *Info, char *name, int *ierr, int namelen )
{
  int i;

  /* printf("Set the controls (for the heart of the sun/pelvis) ...\n"); */

  /* Adjust to 0-based indexing */
  for (i = 0; i <=  *n; ++i) --isr[i];
  for (i = 0; i < *nza; ++i) --isc[i];

  *ierr = umfpack_di_numeric(isr, isc, A, *Symbolic, Numeric, Control, Info);
  if (*ierr != UMFPACK_OK) {
    copy_out(name, "umfpack_numeric", namelen);
    return;
  }

  /* Restore to 1-based indexing */
  for (i = 0; i <=  *n; ++i) ++isr[i];
  for (i = 0; i < *nza; ++i) ++isc[i];

  return;
}


void Umfpack4_Solve( double *A, int *n, int *nza, int *isc, int *isr,
                     double *b, double *x, void **Numeric, double *Control,
                     double *Info, int *ierr )
{
  int i;

  /* Adjust to 0-based indexing */
  for (i = 0; i < *nza; ++i) --isc[i];
  for (i = 0; i <=  *n; ++i) --isr[i];

  /* Solve */
  *ierr = umfpack_di_solve( UMFPACK_At, isr, isc, A, x, b, *Numeric,
                            Control, Info);

  /* Restore to 1-based indexing */
  for (i = 0; i < *nza; ++i) ++isc[i];
  for (i = 0; i <=  *n; ++i) ++isr[i];

  return;
}


void Umfpack4_Free( void **Symbolic, int *free_symbolic,
                    void **Numeric,  int *free_numeric )
{
  if (*free_numeric) {
    umfpack_di_free_numeric(Numeric);
    *Numeric = NULL;
  }

  if (*free_symbolic) {
    umfpack_di_free_symbolic(Symbolic);
    *Symbolic = NULL;
  }

  return;
}
