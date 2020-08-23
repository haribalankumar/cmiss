/* ========================================================================== */
/* === UMFPACK_report_control =============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Prints the control settings.  See umfpack_report_control.h
    for details.
*/

#include "umf_internal.h"

GLOBAL void UMFPACK_report_control
(
    const double Control [UMFPACK_CONTROL]
)
{
    Int prl, nb, irstep ;
    double drow, dcol, relpt, relax, relax2, relax3, alloc_init ;

    if (!Control)
    {
	/* default is to print nothing */
	return ;
    }
    else if (SCALAR_IS_NAN (Control [UMFPACK_PRL]))
    {
	/* default is to print nothing */
	return ;
    }
    else
    {
	prl = (Int) Control [UMFPACK_PRL] ;
    }

    if (prl < 2)
    {
	return ;
    }

    PRINTF (("\n%s, Control:\n", UMFPACK_VERSION)) ;

    /* ---------------------------------------------------------------------- */
    /* run-time options */
    /* ---------------------------------------------------------------------- */

    /* This is a "run-time" option because all four umfpack_* versions */
    /* compiled into the UMFPACK library. */

#ifdef DINT
    PRINTF (("    Matrix entry defined as: double\n")) ;
    PRINTF (("    Int (generic integer) defined as: int\n")) ;
#endif
#ifdef DLONG
    PRINTF (("    Matrix entry defined as: double\n")) ;
    PRINTF (("    Int (generic integer) defined as: long\n")) ;
#endif
#ifdef ZINT
    PRINTF (("    Matrix entry defined as: double complex\n")) ;
    PRINTF (("    Int (generic integer) defined as: int\n")) ;
#endif
#ifdef ZLONG
    PRINTF (("    Matrix entry defined as: double complex\n")) ;
    PRINTF (("    Int (generic integer) defined as: long\n")) ;
#endif

    PRINTF (("    "ID": print level: "ID"\n", (Int) INDEX (UMFPACK_PRL), prl)) ;

    if (SCALAR_IS_NAN (Control [UMFPACK_DENSE_ROW]))
    {
	drow = UMFPACK_DEFAULT_DENSE_ROW ;
    }
    else
    {
	drow = Control [UMFPACK_DENSE_ROW] ;
    }
    PRINTF (("    "ID": dense row parameter:    %g\n",
	(Int) INDEX (UMFPACK_DENSE_ROW), drow)) ;
    PRINTF ((
    "       (\"dense\" rows have    > max (16, (%g)*16*sqrt(n_col)) entries)\n",
	drow)) ;

    if (SCALAR_IS_NAN (Control [UMFPACK_DENSE_COL]))
    {
	dcol = UMFPACK_DEFAULT_DENSE_COL ;
    }
    else
    {
	dcol = Control [UMFPACK_DENSE_COL] ;
    }
    PRINTF (("    "ID": dense column parameter: %g\n",
	(Int) INDEX (UMFPACK_DENSE_COL), dcol)) ;
    PRINTF ((
    "       (\"dense\" columns have > max (16, (%g)*16*sqrt(n_row)) entries)\n",
	dcol)) ;

    if (SCALAR_IS_NAN (Control [UMFPACK_PIVOT_TOLERANCE]))
    {
	relpt = UMFPACK_DEFAULT_PIVOT_TOLERANCE ;
    }
    else
    {
	relpt = Control [UMFPACK_PIVOT_TOLERANCE] ;
    }
    PRINTF (("    "ID": pivot tolerance: %g\n",
	(Int) INDEX (UMFPACK_PIVOT_TOLERANCE), relpt)) ;

    if (SCALAR_IS_NAN (Control [UMFPACK_BLOCK_SIZE]))
    {
	nb = UMFPACK_DEFAULT_BLOCK_SIZE ;
    }
    else
    {
	nb = (Int) Control [UMFPACK_BLOCK_SIZE] ;
    }
    PRINTF (("    "ID": block size for dense matrix kernels: "ID"\n",
	(Int) INDEX (UMFPACK_BLOCK_SIZE), nb)) ;

    if (SCALAR_IS_NAN (Control [UMFPACK_ALLOC_INIT]))
    {
	alloc_init = UMFPACK_DEFAULT_ALLOC_INIT ;
    }
    else
    {
	alloc_init = Control [UMFPACK_ALLOC_INIT] ;
    }
    PRINTF (("    "ID": initial allocation ratio: %g\n",
	(Int) INDEX (UMFPACK_ALLOC_INIT), alloc_init)) ;

    if (SCALAR_IS_NAN (Control [UMFPACK_RELAXED_AMALGAMATION]))
    {
	relax = UMFPACK_DEFAULT_RELAXED_AMALGAMATION ;
    }
    else
    {
	relax = Control [UMFPACK_RELAXED_AMALGAMATION] ;
    }
    PRINTF (("    "ID": relaxed amalgamation parameter:  %g\n",
	(Int) INDEX (UMFPACK_RELAXED_AMALGAMATION), relax)) ;

    if (SCALAR_IS_NAN (Control [UMFPACK_RELAXED2_AMALGAMATION]))
    {
	relax2 = UMFPACK_DEFAULT_RELAXED2_AMALGAMATION ;
    }
    else
    {
	relax2 = Control [UMFPACK_RELAXED2_AMALGAMATION] ;
    }
    PRINTF (("    "ID": relaxed2 amalgamation parameter: %g\n",
	(Int) INDEX (UMFPACK_RELAXED2_AMALGAMATION), relax2)) ;

    if (SCALAR_IS_NAN (Control [UMFPACK_RELAXED3_AMALGAMATION]))
    {
	relax3 = UMFPACK_DEFAULT_RELAXED3_AMALGAMATION ;
    }
    else
    {
	relax3 = Control [UMFPACK_RELAXED3_AMALGAMATION] ;
    }
    PRINTF (("    "ID": relaxed3 amalgamation parameter: %g\n",
	(Int) INDEX (UMFPACK_RELAXED3_AMALGAMATION), relax3)) ;

    if (SCALAR_IS_NAN (Control [UMFPACK_IRSTEP]))
    {
	irstep = UMFPACK_DEFAULT_IRSTEP ;
    }
    else
    {
	irstep = (Int) Control [UMFPACK_IRSTEP] ;
    }
    PRINTF (("    "ID": max iterative refinement steps: "ID"\n\n",
	(Int) INDEX (UMFPACK_IRSTEP), irstep)) ;

    /* ---------------------------------------------------------------------- */
    /* compile-time options */
    /* ---------------------------------------------------------------------- */

    PRINTF ((
	"    The following options can only be changed at compile-time:\n")) ;

    PRINTF (("    "ID": BLAS library used:  ",
	(Int) INDEX (UMFPACK_COMPILED_WITH_BLAS))) ;

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


#ifdef MATLAB_MEX_FILE
    PRINTF (("    "ID": compiled for MATLAB"
    " (uses mxMalloc, mxFree, mxRealloc, and mexPrintf)\n",
	(Int) INDEX (UMFPACK_COMPILED_FOR_MATLAB))) ;
#else
#ifdef MATHWORKS
    PRINTF (("    "ID": compiled for MATLAB, using internal utility routines\n"
    "    (uses utMalloc, utFree, utRealloc, and utPrintf)\n",
	(Int) INDEX (UMFPACK_COMPILED_FOR_MATLAB))) ;
    PRINTF (("    (complex version uses utDivideComplex, utFdlibm_hypot)\n")) ;
#else
    PRINTF (("    "ID": compiled for ANSI C"
    " (uses malloc, free, realloc, and printf)\n",
	(Int) INDEX (UMFPACK_COMPILED_FOR_MATLAB))) ;
#endif
#endif

#ifdef GETRUSAGE
    PRINTF (("    "ID": CPU timer is getrusage.\n",
	(Int) INDEX (UMFPACK_COMPILED_WITH_GETRUSAGE))) ;
#else
    PRINTF (("    "ID": CPU timer is ANSI C clock (may wrap around).\n",
	(Int) INDEX (UMFPACK_COMPILED_WITH_GETRUSAGE))) ;
#endif

#ifndef NDEBUG
    PRINTF (("    "ID": compiled with debugging enabled\n"
    "        This will be exceedingly slow!  ",
	(Int) INDEX (UMFPACK_COMPILED_IN_DEBUG_MODE))) ;
#ifdef MATLAB_MEX_FILE
    PRINTF (("Uses mxAssert.\n")) ;
#else
#ifdef MATHWORKS
    PRINTF (("Uses utAssert.\n")) ;
#else
    PRINTF (("Uses ANSI C assert.\n")) ;
#endif
#endif
#else
    PRINTF (("    "ID": compiled for normal operation (debugging disabled)\n",
	(Int) INDEX (UMFPACK_COMPILED_IN_DEBUG_MODE))) ;
#endif

    PRINTF (("    computer/operating system: %s\n", UMFPACK_ARCHITECTURE)) ;
    PRINTF (("    size of int: %g long: %g Int: %g pointer: %g"
	" double: %g Entry: %g (in bytes)\n\n", (double) sizeof (int),
	(double) sizeof (long), (double) sizeof (Int),
	(double) sizeof (void *), (double) sizeof (double),
	(double) sizeof (Entry))) ;

}

