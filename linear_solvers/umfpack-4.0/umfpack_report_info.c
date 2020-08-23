/* ========================================================================== */
/* === UMFPACK_report_info ================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Prints the Info array.  See umfpack_report_info.h for
    details.
*/

#include "umf_internal.h"

#define PRINT_INFO(format,x) { if ((x) >= 0) PRINTF ((format, x)) ; }

/* RATIO macro uses a double relop, but ignore NaN case: */
#define RATIO(a,b,c) (((b) == 0) ? (c) : (((double) a)/((double) b)))

/* ========================================================================== */
/* === print_ratio ========================================================== */
/* ========================================================================== */

PRIVATE void print_ratio
(
    char *what,
    double estimate,
    double actual
)
{
    if (estimate < 0 && actual < 0)	/* double relop, but ignore Nan case */
    {
	return ;
    }
    PRINTF (("    %-27s", what)) ;
    if (estimate >= 0)			/* double relop, but ignore Nan case */
    {
	PRINTF ((" %20.0f", estimate)) ;
    }
    else
    {
	PRINTF (("                    -")) ;
    }
    if (actual >= 0)			/* double relop, but ignore Nan case */
    {
	PRINTF ((" %20.0f", actual)) ;
    }
    else
    {
	PRINTF (("                    -")) ;
    }
    if (estimate >= 0 && actual >= 0)	/* double relop, but ignore Nan case */
    {
	PRINTF ((" %5.0f%%\n", 100 * RATIO (actual, estimate, 1))) ;
    }
    else
    {
	PRINTF (("      -\n")) ;
    }
}


/* ========================================================================== */

GLOBAL void UMFPACK_report_info
(
    const double Control [UMFPACK_CONTROL],
    const double Info [UMFPACK_INFO]
)
{

    double lnz_est, unz_est, lunz_est, lnz, unz, lunz ;
    Int n_row, n_col, n_inner, prl ;

    /* ---------------------------------------------------------------------- */
    /* get control settings and status to determine what to print */
    /* ---------------------------------------------------------------------- */

    if (!Control)
    {
	prl = UMFPACK_DEFAULT_PRL ;
    }
    else if (SCALAR_IS_NAN (Control [UMFPACK_PRL]))
    {
	prl = UMFPACK_DEFAULT_PRL ;
    }
    else
    {
	prl = (Int) Control [UMFPACK_PRL] ;
    }

    if (!Info || prl < 2)
    {
	/* no output generated if Info is (double *) NULL */
	/* or if prl is less than 2 */
	return ;
    }

    /* ---------------------------------------------------------------------- */
    /* print umfpack version */
    /* ---------------------------------------------------------------------- */

    PRINTF (("\n%s, Info:\n", UMFPACK_VERSION)) ;

    /* ---------------------------------------------------------------------- */
    /* print run-time options */
    /* ---------------------------------------------------------------------- */

#ifdef DINT
    PRINTF (("    matrix entry defined as: double\n")) ;
    PRINTF (("    Int (generic integer) defined as: int\n")) ;
#endif
#ifdef DLONG
    PRINTF (("    matrix entry defined as: double\n")) ;
    PRINTF (("    Int (generic integer) defined as: long\n")) ;
#endif
#ifdef ZINT
    PRINTF (("    matrix entry defined as: double complex\n")) ;
    PRINTF (("    Int (generic integer) defined as: int\n")) ;
#endif
#ifdef ZLONG
    PRINTF (("    matrix entry defined as: double complex\n")) ;
    PRINTF (("    Int (generic integer) defined as: long\n")) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* print compile-time options */
    /* ---------------------------------------------------------------------- */

    PRINTF (("    BLAS library used: ")) ;

#if defined (USE_NO_BLAS)
    PRINTF (("none.  UMFPACK will be slow.\n")) ;
#elif defined (USE_C_BLAS)
    PRINTF (("C-BLAS.\n")) ;
#elif defined (USE_MATLAB_BLAS)
    PRINTF (("built-in MATLAB BLAS.\n")) ;
#elif defined (USE_SUNPERF_BLAS)
    PRINTF (("Sun Performance Library BLAS.\n")) ;
#elif defined (USE_SCSL_BLAS)
    PRINTF (("SGI SCSL BLAS.\n")) ;
#elif defined (USE_FORTRAN_BLAS)
    PRINTF (("Fortran BLAS.\n")) ;
#endif

    PRINTF (("    MATLAB: ")) ;
#ifdef MATLAB_MEX_FILE
    PRINTF (("yes.\n")) ;
#else
#ifdef MATHWORKS
    PRINTF (("yes (using internal ut* routines).\n")) ;
#else
    PRINTF (("no.\n")) ;
#endif
#endif

    PRINTF (("    CPU timer: ")) ;
#ifdef GETRUSAGE
    PRINTF (("getrusage.\n")) ;
#else
    PRINTF (("ANSI C clock (may wrap-around).\n")) ;
#endif

#ifndef NDEBUG
    PRINTF (("    Debugging enabled (umfpack will be exceedingly slow!)\n")) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* print n and nz */
    /* ---------------------------------------------------------------------- */

    n_row = (Int) Info [UMFPACK_NROW] ;
    n_col = (Int) Info [UMFPACK_NCOL] ;
    n_inner = MIN (n_row, n_col) ;

    PRINT_INFO ("    number of rows in input matrix:   "ID"\n", n_row) ;
    PRINT_INFO ("    number of colums in input matrix: "ID"\n", n_col) ;
    PRINT_INFO ("    entries in input matrix:          "ID"\n",
	(Int) Info [UMFPACK_NZ]) ;
    PRINT_INFO ("    memory usage is reported below in "ID"-byte Units\n",
	(Int) Info [UMFPACK_SIZE_OF_UNIT]) ;

    PRINT_INFO ("    size of int:  "ID" bytes\n",
	(Int) Info [UMFPACK_SIZE_OF_INT]) ;
    PRINT_INFO ("    size of long:  "ID" bytes\n",
	(Int) Info [UMFPACK_SIZE_OF_LONG]) ;
    PRINT_INFO ("    size of pointer:  "ID" bytes\n",
	(Int) Info [UMFPACK_SIZE_OF_POINTER]) ;
    PRINT_INFO ("    size of numerical entry:  "ID" bytes\n",
	(Int) Info [UMFPACK_SIZE_OF_ENTRY]) ;

    /* ---------------------------------------------------------------------- */
    /* estimates vs actual nonzeros in L and U */
    /* ---------------------------------------------------------------------- */

    lnz_est = Info [UMFPACK_LNZ_ESTIMATE] ;
    unz_est = Info [UMFPACK_UNZ_ESTIMATE] ;
    if (lnz_est >= 0 && unz_est >= 0)	/* double relop, but ignore NaN case */
    {
	lunz_est = lnz_est + unz_est - n_inner ;
    }
    else
    {
	lunz_est = EMPTY ;
    }
    lnz = Info [UMFPACK_LNZ] ;
    unz = Info [UMFPACK_UNZ] ;
    if (lnz >= 0 && unz >= 0)		/* double relop, but ignore NaN case */
    {
	lunz = lnz + unz - n_inner ;
    }
    else
    {
	lunz = EMPTY ;
    }

    /* ---------------------------------------------------------------------- */
    /* print header for table of estimates/actual statistics */
    /* ---------------------------------------------------------------------- */

    PRINTF (("    symbolic/numeric factorization:         estimate")) ;
    PRINTF (("               actual      %%\n")) ;

    /* ---------------------------------------------------------------------- */
    /* symbolic factorization */
    /* ---------------------------------------------------------------------- */

    PRINT_INFO (
    "    number of \"dense\" rows                         - %20.0f      -\n",
	Info [UMFPACK_NDENSE_ROW]) ;
    PRINT_INFO (
    "    rows with entries only in \"dense\" columns      - %20.0f      -\n",
	Info [UMFPACK_NEMPTY_ROW]) ;
    PRINT_INFO (
    "    number of \"dense\" columns                      - %20.0f      -\n",
	Info [UMFPACK_NDENSE_COL]) ;
    PRINT_INFO (
    "    columns with entries only in \"dense\" rows      - %20.0f      -\n",
	Info [UMFPACK_NEMPTY_COL]) ;
    PRINT_INFO (
    "    symbolic factorization defragmentations        - %20.0f      -\n",
	Info [UMFPACK_SYMBOLIC_DEFRAG]) ;
    PRINT_INFO (
    "    symbolic memory usage (Units)                  - %20.0f      -\n",
	Info [UMFPACK_SYMBOLIC_PEAK_MEMORY]) ;
    PRINT_INFO (
    "    Symbolic size (Units)                          - %20.0f      -\n",
	Info [UMFPACK_SYMBOLIC_SIZE]) ;
    PRINT_INFO (
    "    symbolic factorization time (sec)              - %20.2f      -\n",
	Info [UMFPACK_SYMBOLIC_TIME]) ;

    /* ---------------------------------------------------------------------- */
    /* estimate/actual in symbolic/numeric factorization */
    /* ---------------------------------------------------------------------- */

    /* double relop, but ignore NaN case: */
    if (Info [UMFPACK_SYMBOLIC_DEFRAG] >= 0	/* UMFPACK_*symbolic called */
    ||  Info [UMFPACK_NUMERIC_DEFRAG] >= 0)	/* UMFPACK_numeric called */
    {
	PRINTF (("    variable-sized part of Numeric object:\n")) ;
    }
    print_ratio ("    initial size (Units)",
	Info [UMFPACK_VARIABLE_INIT_ESTIMATE], Info [UMFPACK_VARIABLE_INIT]) ;
    print_ratio ("    peak size (Units)",
	Info [UMFPACK_VARIABLE_PEAK_ESTIMATE], Info [UMFPACK_VARIABLE_PEAK]) ;
    print_ratio ("    final size (Units)",
	Info [UMFPACK_VARIABLE_FINAL_ESTIMATE], Info [UMFPACK_VARIABLE_FINAL]) ;
    print_ratio ("Numeric final size (Units)",
	Info [UMFPACK_NUMERIC_SIZE_ESTIMATE], Info [UMFPACK_NUMERIC_SIZE]) ;
    print_ratio ("peak memory usage (Units)",
	Info [UMFPACK_PEAK_MEMORY_ESTIMATE], Info [UMFPACK_PEAK_MEMORY]) ;
    print_ratio ("numeric factorization flops",
	Info [UMFPACK_FLOPS_ESTIMATE], Info [UMFPACK_FLOPS]) ;
    print_ratio ("nz in L (incl diagonal)", lnz_est, lnz) ;
    print_ratio ("nz in U (incl diagonal)", unz_est, unz) ;
    print_ratio ("nz in L+U (incl diagonal)", lunz_est, lunz) ;
    print_ratio ("largest front (# entries)",
	Info [UMFPACK_MAX_FRONT_SIZE_ESTIMATE], Info [UMFPACK_MAX_FRONT_SIZE]) ;

    /* ---------------------------------------------------------------------- */
    /* numeric factorization */
    /* ---------------------------------------------------------------------- */

    PRINT_INFO (
    "    nonzeros on diagonal of U                      - %20.0f      -\n",
	Info [UMFPACK_UDIAG_NZ]) ;
    PRINT_INFO (
    "    estimate of reciprocal of condition number     - %20.5e      -\n",
	Info [UMFPACK_RCOND]) ;
    PRINT_INFO (
    "    indices in compressed pattern                  - %20.0f      -\n",
	Info [UMFPACK_COMPRESSED_PATTERN]) ;
    PRINT_INFO (
    "    numerical values stored in Numeric object      - %20.0f      -\n",
	Info [UMFPACK_LU_ENTRIES]) ;
    PRINT_INFO (
    "    numeric factorization defragmentations         - %20.0f      -\n",
	Info [UMFPACK_NUMERIC_DEFRAG]) ;
    PRINT_INFO (
    "    numeric factorization reallocations            - %20.0f      -\n",
	Info [UMFPACK_NUMERIC_REALLOC]) ;
    PRINT_INFO (
    "    costly numeric factorization reallocations     - %20.0f      -\n",
	Info [UMFPACK_NUMERIC_COSTLY_REALLOC]) ;
    PRINT_INFO (
    "    numeric factorization time (sec)               - %20.2f      -\n",
	Info [UMFPACK_NUMERIC_TIME]) ;

    /* ---------------------------------------------------------------------- */
    /* solve */
    /* ---------------------------------------------------------------------- */

    PRINT_INFO (
    "    solve flops                                    - %20.0f      -\n",
	Info [UMFPACK_SOLVE_FLOPS]) ;
    PRINT_INFO (
    "    iterative refinement steps taken               - %20.0f      -\n",
	Info [UMFPACK_IR_TAKEN]) ;
    PRINT_INFO (
    "    iterative refinement steps attempted           - %20.0f      -\n",
	Info [UMFPACK_IR_ATTEMPTED]) ;
    PRINT_INFO (
    "    sparse backward error omega1                   - %20g      -\n",
	Info [UMFPACK_OMEGA1]) ;
    PRINT_INFO (
    "    sparse backward error omega2                   - %20g      -\n",
	Info [UMFPACK_OMEGA2]) ;
    PRINT_INFO (
    "    solve time (sec)                               - %20.2f      -\n",
	Info [UMFPACK_SOLVE_TIME]) ;

    PRINTF (("\n")) ;

}

