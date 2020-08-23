/* ========================================================================== */
/* === umfpack mexFunction ================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    MATLAB interface for umfpack.

    Factor or solve a sparse linear system, returning either the solution
    x to Ax=b or A'x'=b', or the factorization LU=PAQ.  A must be sparse, with
    nonzero dimensions, but it may be complex, singular, and/or rectangular.
    b must be a dense n-by-1 vector (real or complex).

    See umfpack.m and umfpack.h for details.

    Note that this mexFunction accesses only the user-callable UMFPACK routines.
    Thus, is also provides another example of how user C code can access
    UMFPACK.

    If NO_TRANSPOSE_FORWARD_SLASH is not defined at compile time, then the
    forward slash (/) operator acts almost like MATLAB's x = b/A.  It is solved
    by factorizing the array transpose, and then x = (A.'\b.').' is solved.
    Note that MATLAB forms the complex conjugate transpose, but this is not
    necessary.  This is the default behavior, since factorizing A can behave
    perform much differently than factorizing its transpose.

    If NO_TRANSPOSE_FORWARD_SLASH is defined at compile time, then the forward
    slash operator does not act like MATLAB's x=b/A.  It is solved by
    factorizing A, and then solving via the transposed L and U matrices.
    The solution is still x = (A.'\b.').', except that A is factorized instead
    of A.'.
*/


/* ========================================================================== */
/* === Include files ======================================================== */
/* ========================================================================== */

#include "umfpack.h"

#include "mex.h"
#include "matrix.h"
#include <string.h>

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define STRING_MATCH(s1,s2) (strcmp ((s1), (s2)) == 0)
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

static void error
(
    char *s,
    int A_is_complex,
    int nlhs,
    mxArray *plhs [ ],
    double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO],
    int status,
    int do_info
)
{
    int i ;
    double *OutInfo ;
    if (A_is_complex)
    {
	umfpack_zi_report_status (Control, status) ;
	umfpack_zi_report_info (Control, Info) ;
    }
    else
    {
	umfpack_di_report_status (Control, status) ;
	umfpack_di_report_info (Control, Info) ;
    }
    if (do_info > 0)
    {
	/* return Info */
	plhs [do_info] = mxCreateDoubleMatrix (1, UMFPACK_INFO, mxREAL) ;
	OutInfo = mxGetPr (plhs [do_info]) ;
	for (i = 0 ; i < UMFPACK_INFO ; i++)
	{
	    OutInfo [i] = Info [i] ;
	}
    }
    mexErrMsgTxt (s) ;
}


/* ========================================================================== */
/* === umfpack ============================================================== */
/* ========================================================================== */

void mexFunction
(
    int nlhs,			/* number of left-hand sides */
    mxArray *plhs [ ],		/* left-hand side matrices */
    int nrhs,			/* number of right--hand sides */
    const mxArray *prhs [ ]	/* right-hand side matrices */
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    void *Symbolic, *Numeric ;
    int *Lp, *Li, *Up, *Ui, *Ap, *Ai, *P, *Q, do_solve, lnz, unz, nn, i,
	transpose, size, do_info, do_numeric, *Qtree, *Front_npivcol, op, k,
	*Front_parent, *Chain_start, *Chain_maxrows, *Chain_maxcols, nz, status,
	nfronts, nchains, *Ltp, *Ltj, *Qinit, print_level, status2,
	*Ptree, *Front_1strow, *Front_leftmostdesc, n_row, n_col, n_inner, sys,
	ignore1, ignore2, ignore3, A_is_complex, B_is_complex, X_is_complex ;
    double *Lx, *Lz, *Ux, *Uz, *Ax, *Az, *Bx, *Bz, *Xx, *Xz, *User_Control,
	*p, *q, Info [UMFPACK_INFO], Control [UMFPACK_CONTROL], *OutInfo,
	*p1, *p2, *p3, *p4, *Ltx, *Ltz ;
    mxArray *Amatrix, *Bmatrix, *User_Control_matrix, *User_Qinit ;
    char *operator, *operation ;
    mxComplexity Atype, Xtype ;
    char warning [200] ;

#ifndef NO_TRANSPOSE_FORWARD_SLASH
    int *Cp, *Ci ;
    double *Cx, *Cz ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get inputs A, b, and the operation to perform */
    /* ---------------------------------------------------------------------- */

    User_Control_matrix = (mxArray *) NULL ;
    User_Qinit = (mxArray *) NULL ;

    do_info = 0 ;
    do_solve = FALSE ;
    do_numeric = TRUE ;
    transpose = FALSE ;

    if (nrhs >= 2 && mxIsChar (prhs [1]))
    {
	op = 1 ;
    }
    else if (nrhs >= 3 && mxIsChar (prhs [2]))
    {
	op = 2 ;
    }
    else
    {
	/* no operator */
	op = 0 ;
    }

    if (op > 0)
    {
	operator = mxArrayToString (prhs [op]) ;

	if (STRING_MATCH (operator, "\\"))
	{

	    /* -------------------------------------------------------------- */
	    /* matrix left divide, x = A\b */
	    /* -------------------------------------------------------------- */

	    /*
		[x, Info] = umfpack (A, '\', b) ;
		[x, Info] = umfpack (A, '\', b, Control) ;
		[x, Info] = umfpack (A, Qinit, '\', b, Control) ;
		[x, Info] = umfpack (A, Qinit, '\', b) ;
	    */

	    operation = "x = A\\b" ;
	    do_solve = TRUE ;
	    Amatrix = (mxArray *) prhs [0] ;
	    Bmatrix = (mxArray *) prhs [op+1] ;

	    if (nlhs == 2)
	    {
		do_info = 1 ;
	    }
	    if (op == 2)
	    {
		User_Qinit = (mxArray *) prhs [1] ;
	    }
	    if ((op == 1 && nrhs == 4) || (op == 2 && nrhs == 5))
	    {
		User_Control_matrix = (mxArray *) prhs [nrhs-1] ;
	    }
	    if (nrhs < 3 || nrhs > 5 || nlhs > 2)
	    {
		mexErrMsgTxt ("wrong number of arguments") ;
	    }

	}
	else if (STRING_MATCH (operator, "/"))
	{

	    /* -------------------------------------------------------------- */
	    /* matrix right divide, x = b/A */
	    /* -------------------------------------------------------------- */

	    /*
		[x, Info] = umfpack (b, '/', A) ;
		[x, Info] = umfpack (b, '/', A, Control) ;
		[x, Info] = umfpack (b, '/', A, Qinit) ;
		[x, Info] = umfpack (b, '/', A, Qinit, Control) ;
	    */

	    operation = "x = b/A" ;
	    do_solve = TRUE ;
	    transpose = TRUE ;
	    Amatrix = (mxArray *) prhs [2] ;
	    Bmatrix = (mxArray *) prhs [0] ;

	    if (nlhs == 2)
	    {
		do_info = 1 ;
	    }
	    if (nrhs == 5)
	    {
		User_Qinit = (mxArray *) prhs [3] ;
		User_Control_matrix = (mxArray *) prhs [4] ;
	    }
	    else if (nrhs == 4)
	    {
		/* Control is k-by-1 where k > 1, Qinit is 1-by-n */
		if (mxGetM (prhs [3]) == 1)
		{
		    User_Qinit = (mxArray *) prhs [3] ;
		}
		else
		{
		    User_Control_matrix = (mxArray *) prhs [3] ;
		}
	    }
	    else if (nrhs < 3 || nrhs > 5 || nlhs > 2)
	    {
		mexErrMsgTxt ("wrong number of arguments") ;
	    }

	}
	else if (STRING_MATCH (operator, "symbolic"))
	{

	    /* -------------------------------------------------------------- */
	    /* symbolic factorization only */
	    /* -------------------------------------------------------------- */

	    /*
	    [Ptree, Qtree, Fr, Ch, Info] = umfpack (A, 'symbolic') ;
	    [Ptree, Qtree, Fr, Ch, Info] = umfpack (A, 'symbolic', Control) ;
	    [Ptree, Qtree, Fr, Ch, Info] = umfpack (A, Qinit, 'symbolic') ;
	    [Ptree, Qtree, Fr, Ch, Info] = umfpack (A, Qinit, 'symbolic',
	    	Control) ;
	    */

	    operation = "symbolic factorization" ;
	    do_numeric = FALSE ;
	    Amatrix = (mxArray *) prhs [0] ;

	    if (nlhs == 5)
	    {
		do_info = 4 ;
	    }
	    if (op == 2)
	    {
		User_Qinit = (mxArray *) prhs [1] ;
	    }
	    if ((op == 1 && nrhs == 3) || (op == 2 && nrhs == 4))
	    {
		User_Control_matrix = (mxArray *) prhs [nrhs-1] ;
	    }
	    if (nrhs < 2 || nrhs > 4 || nlhs > 5 || nlhs < 4)
	    {
		mexErrMsgTxt ("wrong number of arguments") ;
	    }

	}
	else
	{
	    mexErrMsgTxt ("operator must be '/', '\\', or 'symbolic'") ;
	}
	mxFree (operator) ;

    }
    else if (nrhs > 0)
    {

	/* ------------------------------------------------------------------ */
	/* LU factorization */
	/* ------------------------------------------------------------------ */

	/*
	    [L, U, P, Q, Info] = umfpack (A) ;
	    [L, U, P, Q, Info] = umfpack (A, Control) ;
	    [L, U, P, Q, Info] = umfpack (A, Qinit) ;
	    [L, U, P, Q, Info] = umfpack (A, Qinit, Control) ;
	*/

	operation = "numeric factorization" ;
	Amatrix = (mxArray *) prhs [0] ;

	if (nlhs == 5)
	{
	    do_info = 4 ;
	}
	if (nrhs == 3)
	{
	    User_Qinit = (mxArray *) prhs [1] ;
	    User_Control_matrix = (mxArray *) prhs [2] ;
	}
	else if (nrhs == 2)
	{
	    /* Control is k-by-1 where k > 1, Qinit is 1-by-n */
	    if (mxGetM (prhs [1]) == 1)
	    {
		User_Qinit = (mxArray *) prhs [1] ;
	    }
	    else
	    {
		User_Control_matrix = (mxArray *) prhs [1] ;
	    }
	}
	else if (nrhs > 3 || nlhs > 5 || nlhs < 4)
	{
	    mexErrMsgTxt ("wrong number of arguments") ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* return default control settings */
	/* ------------------------------------------------------------------ */

	/*
	    Control = umfpack ;
	    umfpack ;
	*/

	if (nlhs > 1)
	{
	    mexErrMsgTxt ("wrong number of arguments") ;
	}

	plhs [0] = mxCreateDoubleMatrix (UMFPACK_CONTROL, 1, mxREAL) ;
	User_Control = mxGetPr (plhs [0]) ;
	umfpack_di_defaults (User_Control) ;
	if (nlhs == 0)
	{
	    umfpack_di_report_control (User_Control) ;
	}
	return ;
    }

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (mxGetNumberOfDimensions (Amatrix) != 2)
    {
	mexErrMsgTxt ("input matrix A must be 2-dimensional") ;
    }
    n_row = mxGetM (Amatrix) ;
    n_col = mxGetN (Amatrix) ;
    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;
    if (do_solve && n_row != n_col)
    {
	mexErrMsgTxt ("input matrix A must square for '\\' or '/'") ;
    }
    if (!mxIsSparse (Amatrix))
    {
	mexErrMsgTxt ("input matrix A must be sparse") ;
    }

    /* The real/complex status of A determines which version to use, */
    /* (umfpack_di_* or umfpack_zi_*). */
    A_is_complex = mxIsComplex (Amatrix) ;
    Atype = A_is_complex ? mxCOMPLEX : mxREAL ;
    Ap = mxGetJc (Amatrix) ;
    Ai = mxGetIr (Amatrix) ;
    Ax = mxGetPr (Amatrix) ;
    Az = mxGetPi (Amatrix) ;

    if (do_solve)
    {
	if (n_row != n_col)
	{
	    mexErrMsgTxt ("A must be square for \\ or /") ;
	}
	if (transpose)
	{
	    if (mxGetM (Bmatrix) != 1 && mxGetN (Bmatrix) != nn)
	    {
		mexErrMsgTxt ("b has the wrong dimensions") ;
	    }
	}
	else
	{
	    if (mxGetM (Bmatrix) != nn && mxGetN (Bmatrix) != 1)
	    {
		mexErrMsgTxt ("b has the wrong dimensions") ;
	    }
	}
	if (mxGetNumberOfDimensions (Bmatrix) != 2)
	{
	    mexErrMsgTxt ("input matrix b must be 2-dimensional") ;
	}
	if (mxIsSparse (Bmatrix))
	{
	    mexErrMsgTxt ("input matrix b cannot be sparse") ;
	}
	if (mxGetClassID (Bmatrix) != mxDOUBLE_CLASS)
	{
	    mexErrMsgTxt ("input matrix b must be a double precision matrix") ;
	}

	B_is_complex = mxIsComplex (Bmatrix) ;
	Bx = mxGetPr (Bmatrix) ;
	Bz = mxGetPi (Bmatrix) ;

	X_is_complex = A_is_complex || B_is_complex ;
	Xtype = X_is_complex ? mxCOMPLEX : mxREAL ;
    }

    /* ---------------------------------------------------------------------- */
    /* set the Control parameters */
    /* ---------------------------------------------------------------------- */

    if (A_is_complex)
    {
	umfpack_zi_defaults (Control) ;
    }
    else
    {
	umfpack_di_defaults (Control) ;
    }

    if (User_Control_matrix)
    {
	if (mxGetClassID (User_Control_matrix) != mxDOUBLE_CLASS ||
	    mxIsSparse (User_Control_matrix))
	{
	    mexErrMsgTxt ("Control must be a dense real matrix") ;
	}
	size = UMFPACK_CONTROL ;
	size = MIN (size, mxGetNumberOfElements (User_Control_matrix)) ;
	User_Control = mxGetPr (User_Control_matrix) ;
	for (i = 0 ; i < size ; i++)
	{
	    Control [i] = User_Control [i] ;
	}
    }

    if (mxIsNaN (Control [UMFPACK_PRL]))
    {
	print_level = UMFPACK_DEFAULT_PRL ;
    }
    else
    {
	print_level = (int) Control [UMFPACK_PRL] ;
    }

    /* ---------------------------------------------------------------------- */
    /* get Qinit, if present */
    /* ---------------------------------------------------------------------- */

    if (User_Qinit)
    {
	if (mxGetM (User_Qinit) != 1 || mxGetN (User_Qinit) != n_col)
	{
	    mexErrMsgTxt ("Qinit must be 1-by-n_col") ;
	}
	if (mxGetNumberOfDimensions (User_Qinit) != 2)
	{
	    mexErrMsgTxt ("input Qinit must be 2-dimensional") ;
	}
	if (mxIsComplex (User_Qinit))
	{
	    mexErrMsgTxt ("input Qinit must not be complex") ;
	}
	if (mxGetClassID (User_Qinit) != mxDOUBLE_CLASS)
	{
	    mexErrMsgTxt ("input Qinit must be a double matrix") ;
	}
	if (mxIsSparse (User_Qinit))
	{
	    mexErrMsgTxt ("input Qinit must be dense") ;
	}
	Qinit = (int *) mxMalloc (n_col * sizeof (int)) ;
	p = mxGetPr (User_Qinit) ;
	for (k = 0 ; k < n_col ; k++)
	{
	    /* convert from 1-based to 0-based indexing */
	    Qinit [k] = ((int) (p [k])) - 1 ;
	}

    }
    else
    {
	/* umfpack_*_qsymbolic will call colamd to get Qinit. This is the */
	/* same as calling umfpack_*_symbolic with Qinit set to NULL*/
	Qinit = (int *) NULL ;
    }

    /* ---------------------------------------------------------------------- */
    /* report the inputs A and Qinit */
    /* ---------------------------------------------------------------------- */

    if (print_level >= 2)
    {
	/* print the operation */
	mexPrintf ("\numfpack: %s\n", operation) ;
    }

    if (A_is_complex)
    {
	umfpack_zi_report_control (Control) ;
	if (print_level >= 3) mexPrintf ("\nA: ") ;
	(void) umfpack_zi_report_matrix (n_row, n_col, Ap, Ai, Ax, Az,
	    1, Control) ;
	if (Qinit)
	{
	    if (print_level >= 3) mexPrintf ("\nQinit: ") ;
	    (void) umfpack_zi_report_perm (n_col, Qinit, Control) ;
	}
    }
    else
    {
	umfpack_di_report_control (Control) ;
	if (print_level >= 3) mexPrintf ("\nA: ") ;
	(void) umfpack_di_report_matrix (n_row, n_col, Ap, Ai, Ax,
	    1, Control) ;
	if (Qinit)
	{
	    if (print_level >= 3) mexPrintf ("\nQinit: ") ;
	    (void) umfpack_di_report_perm (n_col, Qinit, Control) ;
	}
    }

#ifndef NO_TRANSPOSE_FORWARD_SLASH
    /* ---------------------------------------------------------------------- */
    /* create the array transpose for x = b/A */
    /* ---------------------------------------------------------------------- */

    if (transpose)
    {
	/* note that in this case A will be square (nn = n_row = n_col) */
	/* x = (A.'\b.').' will be computed */

	/* make sure Ci and Cx exist, avoid malloc of zero-sized arrays. */
	nz = MAX (Ap [nn], 1) ;

	Cp = (int *) mxMalloc ((nn+1) * sizeof (int)) ;
	Ci = (int *) mxMalloc (nz * sizeof (int)) ;
	Cx = (double *) mxMalloc (nz * sizeof (double)) ;
	if (A_is_complex)
	{
	    Cz = (double *) mxMalloc (nz * sizeof (double)) ;
	    status = umfpack_zi_transpose (nn, nn, Ap, Ai, Ax, Az,
	        (int *) NULL, (int *) NULL, Cp, Ci, Cx, Cz, FALSE) ;
	}
	else
	{
	    status = umfpack_di_transpose (nn, nn, Ap, Ai, Ax,
	        (int *) NULL, (int *) NULL, Cp, Ci, Cx) ;
	}

	if (status != UMFPACK_OK)
	{
	    error ("transpose of A failed",
		A_is_complex, nlhs, plhs, Control, Info, status, do_info) ;
	    return ;
	}

	/* modify pointers so that C will be factorized and solved, not A */
	Ap = Cp ;
	Ai = Ci ;
	Ax = Cx ;
	if (A_is_complex)
	{
	    Az = Cz ;
	}
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* perform the symbolic factorization */
    /* ---------------------------------------------------------------------- */

    if (A_is_complex)
    {
	status = umfpack_zi_qsymbolic (n_row, n_col, Ap, Ai, Qinit, &Symbolic,
	    Control, Info) ;
    }
    else
    {
	status = umfpack_di_qsymbolic (n_row, n_col, Ap, Ai, Qinit, &Symbolic,
	    Control, Info) ;
    }

    if (Qinit)
    {
	mxFree (Qinit) ;
    }

    if (status < 0)
    {
	error ("symbolic factorization failed",
	    A_is_complex, nlhs, plhs, Control, Info, status, do_info) ;
	return ;
    }

    /* ---------------------------------------------------------------------- */
    /* report the Symbolic object */
    /* ---------------------------------------------------------------------- */

    if (A_is_complex)
    {
	(void) umfpack_zi_report_symbolic (Symbolic, Control) ;
    }
    else
    {
	(void) umfpack_di_report_symbolic (Symbolic, Control) ;
    }

    /* ---------------------------------------------------------------------- */
    /* perform numeric factorization, or just return symbolic factorization */
    /* ---------------------------------------------------------------------- */

    if (do_numeric)
    {

	/* ------------------------------------------------------------------ */
	/* perform the numeric factorization */
	/* ------------------------------------------------------------------ */

	if (A_is_complex)
	{
	    status = umfpack_zi_numeric (Ap, Ai, Ax, Az, Symbolic, &Numeric,
		Control, Info) ;
	}
	else
	{
	    status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric,
		Control, Info) ;
	}

	/* ------------------------------------------------------------------ */
	/* free the symbolic factorization */
	/* ------------------------------------------------------------------ */

	if (A_is_complex)
	{
	    umfpack_zi_free_symbolic (&Symbolic) ;
	}
	else
	{
	    umfpack_di_free_symbolic (&Symbolic) ;
	}

	/* ------------------------------------------------------------------ */
	/* report the Numeric object */
	/* ------------------------------------------------------------------ */

	if (status < 0)
	{
	    error ("numeric factorization failed",
		A_is_complex, nlhs, plhs, Control, Info, status, do_info) ;
	    return ;
	}

	if (A_is_complex)
	{
	    (void) umfpack_zi_report_numeric (Numeric, Control) ;
	}
	else
	{
	    (void) umfpack_di_report_numeric (Numeric, Control) ;
	}

	/* ------------------------------------------------------------------ */
	/* return the solution or the factorization */
	/* ------------------------------------------------------------------ */

	if (do_solve)
	{
	    /* -------------------------------------------------------------- */
	    /* solve Ax=b or A'x'=b', and return just the solution x */
	    /* -------------------------------------------------------------- */

#ifndef NO_TRANSPOSE_FORWARD_SLASH
	    if (transpose)
	    {
		/* A.'x.'=b.' gives the same x=b/A as solving A'x'=b' */
		/* since C=A.' was factorized, solve with sys = UMFPACK_A */
		/* since x and b are vectors, x.' and b.' are implicit */
		plhs [0] = mxCreateDoubleMatrix (1, nn, Xtype) ;
	    }
	    else
	    {
		plhs [0] = mxCreateDoubleMatrix (nn, 1, Xtype) ;
	    }
	    sys = UMFPACK_A ;
#else
	    if (transpose)
	    {
		/* If A is real, A'x=b is the same as A.'x=b. */
		/* x and b are vectors, so x and b are the same as x' and b'. */
		/* If A is complex, then A.'x.'=b.' gives the same solution x */
		/* as the complex conjugate transpose.  If we used the A'x=b */
		/* option in umfpack_*_solve, we would have to form b' on */
		/* input and x' on output (negating the imaginary part). */
		/* We can save this work by just using the A.'x=b option in */
		/* umfpack_*_solve.  Then, forming x.' and b.' is implicit, */
		/* since x and b are just vectors anyway. */
		/* In both cases, the system to solve is A.'x=b */
		plhs [0] = mxCreateDoubleMatrix (1, nn, Xtype) ;
		sys = UMFPACK_Aat ;
	    }
	    else
	    {
		plhs [0] = mxCreateDoubleMatrix (nn, 1, Xtype) ;
		sys = UMFPACK_A ;
	    }
#endif

	    /* -------------------------------------------------------------- */
	    /* print the right-hand-side, B */
	    /* -------------------------------------------------------------- */

	    if (print_level >= 3) mexPrintf ("\nright-hand side, b: ") ;
	    if (B_is_complex)
	    {
		(void) umfpack_zi_report_vector (nn, Bx, Bz, Control) ;
	    }
	    else
	    {
		(void) umfpack_di_report_vector (nn, Bx, Control) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* solve the system */
	    /* -------------------------------------------------------------- */

	    Xx = mxGetPr (plhs [0]) ;
	    Xz = mxGetPi (plhs [0]) ;
	    status2 = UMFPACK_OK ;

	    if (A_is_complex)
	    {
		if (!B_is_complex)
		{
		    /* umfpack_zi_solve expects a complex B */
		    Bz = (double *) mxCalloc (nn, sizeof (double)) ;
		}
		status = umfpack_zi_solve (sys, Ap, Ai, Ax, Az, Xx, Xz, Bx, Bz,
		    Numeric, Control, Info) ;
		if (!B_is_complex)
		{
		    mxFree (Bz) ;
		}
	    }
	    else
	    {
		if (B_is_complex)
		{
		    /* Ax=b when b is complex and A is sparse can be split */
		    /* into two systems, A*xr=br and A*xi=bi, where r denotes */
		    /* the real part and i the imaginary part of x and b. */
		    status2 = umfpack_di_solve (sys, Ap, Ai, Ax, Xz, Bz,
		    Numeric, Control, Info) ;
		}
		status = umfpack_di_solve (sys, Ap, Ai, Ax, Xx, Bx,
		    Numeric, Control, Info) ;
	    }

#ifndef NO_TRANSPOSE_FORWARD_SLASH
	    /* -------------------------------------------------------------- */
	    /* free the transposed matrix C */
	    /* -------------------------------------------------------------- */

	    if (transpose)
	    {
	        mxFree (Cp) ;
	        mxFree (Ci) ;
	        mxFree (Cx) ;
	        if (A_is_complex)
	        {
	            mxFree (Cz) ;
	        }
	    }
#endif

	    /* -------------------------------------------------------------- */
	    /* free the Numeric object */
	    /* -------------------------------------------------------------- */

	    if (A_is_complex)
	    {
		umfpack_zi_free_numeric (&Numeric) ;
	    }
	    else
	    {
		umfpack_di_free_numeric (&Numeric) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* check error status */
	    /* -------------------------------------------------------------- */

	    if (status < 0 || status2 < 0)
	    {
		mxDestroyArray (plhs [0]) ;
		error ("solve failed",
		    A_is_complex, nlhs, plhs, Control, Info, status, do_info) ;
		return ;
	    }

	    /* -------------------------------------------------------------- */
	    /* print the solution, X */
	    /* -------------------------------------------------------------- */

	    if (print_level >= 3) mexPrintf ("\nsolution, x: ") ;
	    if (X_is_complex)
	    {
		(void) umfpack_zi_report_vector (nn, Xx, Xz, Control) ;
	    }
	    else
	    {
		(void) umfpack_di_report_vector (nn, Xx, Control) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* warn about singular or near-singular matrices */
	    /* -------------------------------------------------------------- */

	    if (status == UMFPACK_WARNING_singular_matrix)
	    {
	    	;
/*
		mexWarnMsgTxt ("matrix is singular") ;
*/
	    }
	    else if (Info [UMFPACK_RCOND] < 1e-14)
	    {
		sprintf (warning, "matrix is nearly singular, rcond = %g",
		    Info [UMFPACK_RCOND]) ;
		mexWarnMsgTxt (warning) ;
	    }

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* get L, U, P, and Q */
	    /* -------------------------------------------------------------- */

	    if (A_is_complex)
	    {
	        status = umfpack_zi_get_lunz (&lnz, &unz, &ignore1, &ignore2,
		    &ignore3, Numeric) ;
	    }
	    else
	    {
	        status = umfpack_di_get_lunz (&lnz, &unz, &ignore1, &ignore2,
		    &ignore3, Numeric) ;
	    }

	    if (status < 0)
	    {
		if (A_is_complex)
		{
		    umfpack_zi_free_numeric (&Numeric) ;
		}
		else
		{
		    umfpack_di_free_numeric (&Numeric) ;
		}
		error ("extracting LU factors failed",
		    A_is_complex, nlhs, plhs, Control, Info, status, do_info) ;
		return ;
	    }

	    /* avoid malloc of zero-sized arrays */
	    lnz = MAX (lnz, 1) ;
	    unz = MAX (unz, 1) ;

	    /* temporary space, for the *** ROW *** form of L */
	    Ltp = (int *) mxMalloc ((n_row+1) * sizeof (int)) ;
	    Ltj = (int *) mxMalloc (lnz * sizeof (int)) ;
	    Ltx = (double *) mxMalloc (lnz * sizeof (double)) ;

	    if (A_is_complex)
	    {
	        Ltz = (double *) mxMalloc (lnz * sizeof (double)) ;
	    }
	    else
	    {
	        Ltz = (double *) NULL ;
	    }

	    plhs [1] = mxCreateSparse (n_inner, n_col, unz, Atype) ;
	    Up = mxGetJc (plhs [1]) ;
	    Ui = mxGetIr (plhs [1]) ;
	    Ux = mxGetPr (plhs [1]) ;
	    Uz = mxGetPi (plhs [1]) ;

	    P = (int *) mxMalloc (n_row * sizeof (int)) ;
	    Q = (int *) mxMalloc (n_col * sizeof (int)) ;

	    if (A_is_complex)
	    {
		status = umfpack_zi_get_numeric (Ltp, Ltj, Ltx, Ltz, Up, Ui,
		    Ux, Uz, P, Q, (double *) NULL, (double *) NULL, Numeric) ;
		umfpack_zi_free_numeric (&Numeric) ;
	    }
	    else
	    {
		status = umfpack_di_get_numeric (Ltp, Ltj, Ltx, Up, Ui,
		    Ux, P, Q, (double *) NULL, Numeric) ;
		umfpack_di_free_numeric (&Numeric) ;
	    }

	    if (status < 0)
	    {
		mxFree (Ltp) ;
		mxFree (Ltj) ;
		mxFree (Ltx) ;
		if (Ltz) mxFree (Ltz) ;
		mxFree (P) ;
		mxFree (Q) ;
		mxDestroyArray (plhs [1]) ;
		error ("extracting LU factors failed",
		    A_is_complex, nlhs, plhs, Control, Info, status, do_info) ;
		return ;
	    }

	    plhs [2] = mxCreateDoubleMatrix (1, n_row, mxREAL) ;
	    plhs [3] = mxCreateDoubleMatrix (1, n_col, mxREAL) ;

	    p = mxGetPr (plhs [2]) ;
	    for (i = 0 ; i < n_row ; i++)
	    {
		/* UMFPACK is 0-based, but MATLAB expects this to be 1-based */
		p [i] = (double) (P [i] + 1) ;
	    }

	    q = mxGetPr (plhs [3]) ;
	    for (i = 0 ; i < n_col ; i++)
	    {
		/* UMFPACK is 0-based, but MATLAB expects this to be 1-based */
		q [i] = (double) (Q [i] + 1) ;
	    }

	    /* permanent copy of L */
	    plhs [0] = mxCreateSparse (n_row, n_inner, lnz, Atype) ;
	    Lp = mxGetJc (plhs [0]) ;
	    Li = mxGetIr (plhs [0]) ;
	    Lx = mxGetPr (plhs [0]) ;
	    Lz = mxGetPi (plhs [0]) ;

	    /* convert L from row form to column form */
	    if (A_is_complex)
	    {
		/* non-conjugate array transpose */
	        status = umfpack_zi_transpose (n_inner, n_row, Ltp, Ltj, Ltx,
		    Ltz, (int *) NULL, (int *) NULL, Lp, Li, Lx, Lz, FALSE) ;
	    }
	    else
	    {
	        status = umfpack_di_transpose (n_inner, n_row, Ltp, Ltj, Ltx,
		    (int *) NULL, (int *) NULL, Lp, Li, Lx) ;
	    }

	    mxFree (Ltp) ;
	    mxFree (Ltj) ;
	    mxFree (Ltx) ;
	    if (Ltz) mxFree (Ltz) ;

	    if (status < 0)
	    {
		mxFree (P) ;
		mxFree (Q) ;
		mxDestroyArray (plhs [0]) ;
		mxDestroyArray (plhs [1]) ;
		mxDestroyArray (plhs [2]) ;
		mxDestroyArray (plhs [3]) ;
		error ("constructing L failed",
		    A_is_complex, nlhs, plhs, Control, Info, status, do_info) ;
		return ;
	    }

	    /* -------------------------------------------------------------- */
	    /* print L, U, P, and Q */
	    /* -------------------------------------------------------------- */

	    if (A_is_complex)
	    {
		if (print_level >= 3) mexPrintf ("\nL: ") ;
	        (void) umfpack_zi_report_matrix (n_row, n_inner, Lp, Li,
		    Lx, Lz, 1, Control) ;
		if (print_level >= 3) mexPrintf ("\nU: ") ;
	        (void) umfpack_zi_report_matrix (n_inner, n_col,  Up, Ui,
		    Ux, Uz, 1, Control) ;
		if (print_level >= 3) mexPrintf ("\nP: ") ;
	        (void) umfpack_zi_report_perm (n_row, P, Control) ;
		if (print_level >= 3) mexPrintf ("\nQ: ") ;
	        (void) umfpack_zi_report_perm (n_col, Q, Control) ;
	    }
	    else
	    {
		if (print_level >= 3) mexPrintf ("\nL: ") ;
	        (void) umfpack_di_report_matrix (n_row, n_inner, Lp, Li,
		    Lx, 1, Control) ;
		if (print_level >= 3) mexPrintf ("\nU: ") ;
	        (void) umfpack_di_report_matrix (n_inner, n_col,  Up, Ui,
		    Ux, 1, Control) ;
		if (print_level >= 3) mexPrintf ("\nP: ") ;
	        (void) umfpack_di_report_perm (n_row, P, Control) ;
		if (print_level >= 3) mexPrintf ("\nQ: ") ;
	        (void) umfpack_di_report_perm (n_col, Q, Control) ;
	    }

	    mxFree (P) ;
	    mxFree (Q) ;

	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* return the symbolic factorization */
	/* ------------------------------------------------------------------ */

	plhs [0] = mxCreateDoubleMatrix (1, n_row, mxREAL) ;
	plhs [1] = mxCreateDoubleMatrix (1, n_col, mxREAL) ;
	Qtree = (int *) mxMalloc (n_col * sizeof (int)) ;
	Ptree = (int *) mxMalloc (n_row * sizeof (int)) ;
	Front_npivcol = (int *) mxMalloc ((nn+1) * sizeof (int)) ;
	Front_parent = (int *) mxMalloc ((nn+1) * sizeof (int)) ;
	Front_1strow = (int *) mxMalloc ((nn+1) * sizeof (int)) ;
	Front_leftmostdesc = (int *) mxMalloc ((nn+1) * sizeof (int)) ;
	Chain_start = (int *) mxMalloc ((nn+1) * sizeof (int)) ;
	Chain_maxrows = (int *) mxMalloc ((nn+1) * sizeof (int)) ;
	Chain_maxcols = (int *) mxMalloc ((nn+1) * sizeof (int)) ;

	if (A_is_complex)
	{
	    status = umfpack_zi_get_symbolic (&ignore1, &ignore2,
	        &nz, &nfronts, &nchains, Ptree, Qtree, Front_npivcol,
	        Front_parent, Front_1strow, Front_leftmostdesc,
	        Chain_start, Chain_maxrows, Chain_maxcols, Symbolic) ;
	    umfpack_zi_free_symbolic (&Symbolic) ;
	}
	else
	{
	    status = umfpack_di_get_symbolic (&ignore1, &ignore2,
	        &nz, &nfronts, &nchains, Ptree, Qtree, Front_npivcol,
	        Front_parent, Front_1strow, Front_leftmostdesc,
	        Chain_start, Chain_maxrows, Chain_maxcols, Symbolic) ;
	    umfpack_di_free_symbolic (&Symbolic) ;
	}

	if (status < 0)
	{
	    mxFree (Ptree) ;
	    mxFree (Qtree) ;
	    mxFree (Front_npivcol) ;
	    mxFree (Front_parent) ;
	    mxFree (Front_1strow) ;
	    mxFree (Front_leftmostdesc) ;
	    mxFree (Chain_start) ;
	    mxFree (Chain_maxrows) ;
	    mxFree (Chain_maxcols) ;
	    mxDestroyArray (plhs [0]) ;
	    mxDestroyArray (plhs [1]) ;
	    error ("extracting symbolic factors failed",
		A_is_complex, nlhs, plhs, Control, Info, status, do_info) ;
	    return ;
	}

	plhs [2] = mxCreateDoubleMatrix (nfronts+1, 4, mxREAL) ;
	plhs [3] = mxCreateDoubleMatrix (nchains+1, 3, mxREAL) ;

	p = mxGetPr (plhs [0]) ;
	for (i = 0 ; i < n_row ; i++)
	{
	    /* UMFPACK is 0-based, but MATLAB expects this to be 1-based */
	    p [i] = (double) (Ptree [i] + 1) ;
	}

	p = mxGetPr (plhs [1]) ;
	for (i = 0 ; i < n_col ; i++)
	{
	    /* UMFPACK is 0-based, but MATLAB expects this to be 1-based */
	    p [i] = (double) (Qtree [i] + 1) ;
	}

	p1 = mxGetPr (plhs [2]) ;
	p2 = p1 + nfronts + 1 ;
	p3 = p2 + nfronts + 1 ;
	p4 = p3 + nfronts + 1 ;
	for (i = 0 ; i <= nfronts ; i++)
	{
	    /* convert parent, 1strow, and leftmostdesc to 1-based */
	    p1 [i] = (double) (Front_npivcol [i]) ;
	    p2 [i] = (double) (Front_parent [i] + 1) ;
	    p3 [i] = (double) (Front_1strow [i] + 1) ;
	    p4 [i] = (double) (Front_leftmostdesc [i] + 1) ;
	}

	p1 = mxGetPr (plhs [3]) ;
	p2 = p1 + nchains + 1 ;
	p3 = p2 + nchains + 1 ;
	for (i = 0 ; i < nchains ; i++)
	{
	    p1 [i] = (double) (Chain_start [i] + 1) ;	/* convert to 1-based */
	    p2 [i] = (double) (Chain_maxrows [i]) ;
	    p3 [i] = (double) (Chain_maxcols [i]) ;
	}
	p1 [nchains] = Chain_start [nchains] + 1 ;
	p2 [nchains] = 0 ;
	p3 [nchains] = 0 ;

	mxFree (Ptree) ;
	mxFree (Qtree) ;
	mxFree (Front_npivcol) ;
	mxFree (Front_parent) ;
	mxFree (Front_1strow) ;
	mxFree (Front_leftmostdesc) ;
	mxFree (Chain_start) ;
	mxFree (Chain_maxrows) ;
	mxFree (Chain_maxcols) ;
    }

    /* ---------------------------------------------------------------------- */
    /* report Info */
    /* ---------------------------------------------------------------------- */

    if (A_is_complex)
    {
	umfpack_zi_report_info (Control, Info) ;
    }
    else
    {
	umfpack_di_report_info (Control, Info) ;
    }

    if (do_info > 0)
    {
	/* return Info */
	plhs [do_info] = mxCreateDoubleMatrix (1, UMFPACK_INFO, mxREAL) ;
	OutInfo = mxGetPr (plhs [do_info]) ;
	for (i = 0 ; i < UMFPACK_INFO ; i++)
	{
	    OutInfo [i] = Info [i] ;
	}
    }
}

