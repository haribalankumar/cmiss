/* ========================================================================== */
/* === UMFPACK_report_numeric =============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Prints the Numeric object.
    See umfpack_report_numeric.h for details.

    Dynamic memory usage:  Allocates a size n*sizeof(Int)
    workspace via a single call to UMF_malloc and then frees all of it via
    UMF_free on return.  The workspace is not allocated if an early error
    return occurs  before the workspace is needed.

*/

#include "umf_internal.h"
#include "umf_valid_numeric.h"
#include "umf_report_perm.h"
#include "umf_report_vector.h"
#include "umf_malloc.h"
#include "umf_free.h"


PRIVATE Int UMF_report_L
(
    NumericType *Numeric,
    Int Pattern [ ],
    Int prl
) ;


PRIVATE Int UMF_report_U
(
    NumericType *Numeric,
    Int Pattern [ ],
    Int prl
) ;


GLOBAL Int UMFPACK_report_numeric
(
    void *NumericHandle,
    const double Control [UMFPACK_CONTROL]
)
{
    Int prl, *W, nn, n_row, n_col, n_inner, num_fixed_size, numeric_size ;
    NumericType *Numeric ;

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

    if (prl <= 2)
    {
	return (UMFPACK_OK) ;
    }

    PRINTF (("Numeric object:  ")) ;

    Numeric = (NumericType *) NumericHandle ;
    if (!UMF_valid_numeric (Numeric))
    {
	PRINTF (("ERROR: invalid\n\n")) ;
	return (UMFPACK_ERROR_invalid_Numeric_object) ;
    }

    n_row = Numeric->n_row ;
    n_col = Numeric->n_col ;
    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;

    /* size of Numeric object, except Numeric->Memory and Numeric->Upattern */
    /* see also UMF_set_stats */
    num_fixed_size =
	UNITS (NumericType, 1)		/* Numeric structure */
	+ UNITS (Entry, n_inner+1)	/* D */
	+ 4 * UNITS (Int, n_row+1)	/* Rperm, Lpos, Uilen, Uip */
	+ 4 * UNITS (Int, n_col+1) ;	/* Cperm, Upos, Lilen, Lip */
    
    /* size of Numeric->Memory is Numeric->size */
    /* size of Numeric->Upattern is Numeric->ulen */
    numeric_size = num_fixed_size + Numeric->size + Numeric->ulen ;

    if (prl >= 4)
    {
	PRINTF (("\n    n_row: "ID"  n_col: "ID"\n", n_row, n_col)) ;
	PRINTF (("    relative pivot tolerance used:              %g\n",
	    Numeric->relpt)) ;
	PRINTF (("    relaxed amalgamation parameters used:       %g %g %g\n",
	    Numeric->relax, Numeric->relax2, Numeric->relax3)) ;
	PRINTF (("    initial allocation parameter used:          %g\n",
	    Numeric->alloc_init)) ;
	PRINTF (("    final total size of Numeric object (Units): "ID"\n",
	    numeric_size)) ;
	PRINTF (("    peak size of variable-size part (Units):    "ID"\n",
	    Numeric->max_usage)) ;
	PRINTF (("    largest actual frontal matrix size:         "ID"\n",
	    Numeric->maxfrsize)) ;
	PRINTF (("    memory defragmentations:                    "ID"\n",
	    Numeric->ngarbage)) ;
	PRINTF (("    memory reallocations:                       "ID"\n",
	    Numeric->nrealloc)) ;
	PRINTF (("    costly memory reallocations:                "ID"\n",
	    Numeric->ncostly)) ;
	PRINTF (("    entries in compressed pattern (L and U):    "ID"\n",
	    Numeric->isize)) ;
	PRINTF (("    number of nonzeros in L (excl diag):        "ID"\n",
	    Numeric->lnz)) ;
	PRINTF (("    number of entries stored in L (excl diag):  "ID"\n",
	    Numeric->nLentries)) ;
	PRINTF (("    number of nonzeros in U (excl diag):        "ID"\n",
	    Numeric->unz)) ;
	PRINTF (("    number of entries stored in U (excl diag):  "ID"\n",
	    Numeric->nUentries)) ;
	PRINTF (("    factorization floating-point operations:    %g\n",
	    Numeric->flops)) ;
	PRINTF (("    number of nonzeros on diagonal of U:        "ID"\n",
	    Numeric->nnzpiv)) ;
	PRINTF (("    smallest abs. value on diagonal of U:       %g\n",
	    Numeric->min_udiag)) ;
	PRINTF (("    largest abs. value on diagonal of U:        %g\n",
	    Numeric->max_udiag)) ;
	PRINTF (("    reciprical condition number estimate:       %g\n",
	    Numeric->rcond)) ;
    }

    W = (Int *) UMF_malloc (nn, sizeof (Int)) ;
    if (!W)
    {
	PRINTF ((" ERROR: out of memory to check Numeric object\n\n")) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }

    PRINTF4 (("\nP: row ")) ;
    if (UMF_report_perm (n_row, Numeric->Rperm, W, prl, 0) != UMFPACK_OK)
    {
	(void) UMF_free ((void *) W) ;
	return (UMFPACK_ERROR_invalid_Numeric_object) ;
    }

    PRINTF4 (("\nQ: column ")) ;
    if (UMF_report_perm (n_col, Numeric->Cperm, W, prl, 0) != UMFPACK_OK)
    {
	(void) UMF_free ((void *) W) ;
	return (UMFPACK_ERROR_invalid_Numeric_object) ;
    }

    if (!UMF_report_L (Numeric, W, prl))
    {
	(void) UMF_free ((void *) W) ;
	PRINTF ((" ERROR: invalid\n\n")) ;
	return (UMFPACK_ERROR_invalid_Numeric_object) ;
    }

    if (!UMF_report_U (Numeric, W, prl))
    {
	(void) UMF_free ((void *) W) ;
	PRINTF ((" ERROR: invalid\n\n")) ;
	return (UMFPACK_ERROR_invalid_Numeric_object) ;
    }

    /* The diagonal of U is in "merged" (Entry) form, not "split" form. */
    PRINTF4 (("\nDiagonal of U: ")) ;
    (void) UMF_report_vector (n_inner, (double *) Numeric->D, (double *) NULL,
	prl, 0) ;

    (void) UMF_free ((void *) W) ;

    PRINTF4 (("    Numeric object:  ")) ;
    PRINTF (("OK\n\n")) ;
    return (UMFPACK_OK) ;
}


/* ========================================================================== */
/* === UMF_report_L ========================================================= */
/* ========================================================================== */

PRIVATE Int UMF_report_L
(
    NumericType *Numeric,
    Int Pattern [ ],
    Int prl
)
{
    Int k, deg, *ip, j, row, n_row, *Lpos, *Lilen, valid, k1,
	*Lip, newLchain, llen, prl1, pos, lp, p, npiv ;
    Entry *xp ;

    /* ---------------------------------------------------------------------- */

    ASSERT (prl >= 3) ;

    n_row = Numeric->n_row ;
    npiv = Numeric->npiv ;
    Lpos = Numeric->Lpos ;
    Lilen = Numeric->Lilen ;
    Lip = Numeric->Lip ;
    prl1 = prl ;
    deg = 0 ;

    PRINTF4 ((
    "\nL in Numeric object, in column-oriented compressed-pattern form:\n"
    "    Diagonal entries are all equal to 1.0 (not stored)\n")) ;

    ASSERT (Pattern) ;

    /* ---------------------------------------------------------------------- */
    /* print L */
    /* ---------------------------------------------------------------------- */

    k1 = 12 ;

    for (k = 0 ; k < npiv ; k++)
    {
	/* if prl is 4, print the first 10 entries of the first 10 columns */
	if (k1 > 0)
	{
	    prl = prl1 ;
	}

	lp = Lip [k] ;
	newLchain = (lp < 0) ;
	if (newLchain)
	{
	    lp = -lp ;
	    deg = 0 ;
	}

	if (k1-- > 0)
	{
	    prl = prl1 ;
	}
	else if (prl == 4)
	{
	    PRINTF (("    ...\n")) ;
	    prl-- ;
	}

	PRINTF4 (("\n    column "ID":", INDEX (k))) ;

	/* ------------------------------------------------------------------ */
	/* make column of L in Pattern [0..deg-1] */
	/* ------------------------------------------------------------------ */

	/* remove pivot row */
	pos = Lpos [k] ;
	if (pos != EMPTY)
	{
	    PRINTF4 (("  remove row "ID" at position "ID".",
		INDEX (Pattern [pos]), INDEX (pos))) ;
	    valid = (!newLchain) && (deg > 0) && (pos < deg) && (pos >= 0)
		&& (Pattern [pos] == k) ;
	    if (!valid)
	    {
		return (FALSE) ;
	    }
	    Pattern [pos] = Pattern [--deg] ;
	}

	/* concatenate the pattern */
	llen = Lilen [k] ;
	if (llen < 0)
	{
	    return (FALSE) ;
	}
	p = lp + UNITS (Int, llen) ;
	xp = (Entry *) (Numeric->Memory + p) ;
	if ((llen > 0 || deg > 0)
	    && (p + (Int) UNITS (Entry, deg) > Numeric->size))
	{
	    return (FALSE) ;
	}
	if (llen > 0)
	{
	    PRINTF4 (("  add "ID" entries.", llen)) ;
	    ip = (Int *) (Numeric->Memory + lp) ;
	    for (j = 0 ; j < llen ; j++)
	    {
		Pattern [deg++] = *ip++ ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* print column k of L */
	/* ------------------------------------------------------------------ */

	PRINTF4 (("  length "ID".", deg)) ;
	if (newLchain)
	{
	    PRINTF4 (("  Start of Lchain.")) ;
	}
	PRINTF4 (("\n")) ;

	for (j = 0 ; j < deg ; j++)
	{
	    row = Pattern [j] ;
	    PRINTF4 (("\trow "ID" : ", INDEX (row))) ;
	    if (prl >= 4) PRINT_ENTRY (*xp) ;
	    if (row <= k || row >= n_row)
	    {
		return (FALSE) ;
	    }
	    PRINTF4 (("\n")) ;
	    xp++ ;
	    /* truncate printout, but continue to check L */
	    if (prl == 4 && j == 9 && deg > 10)
	    {
		PRINTF (("\t...\n")) ;
		prl-- ;
	    }
	}
    }

    PRINTF4 (("\n")) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === UMF_report_U ========================================================= */
/* ========================================================================== */

PRIVATE Int UMF_report_U
(
    NumericType *Numeric,
    Int Pattern [ ],
    Int prl
)
{
    /* ---------------------------------------------------------------------- */

    Int k, deg, j, *ip, col, *Upos, *Uilen, k1, prl1, pos,
	*Uip, n_col, ulen, p, newUchain, up, npiv ;
    Entry *xp ;

    /* ---------------------------------------------------------------------- */

    ASSERT (prl >= 3) ;

    n_col = Numeric->n_col ;
    npiv = Numeric->npiv ;
    Upos = Numeric->Upos ;
    Uilen = Numeric->Uilen ;
    Uip = Numeric->Uip ;
    prl1 = prl ;

    PRINTF4 ((
    "\nU in Numeric object, in row-oriented compressed-pattern form:\n"
    "    Diagonal is stored separately.\n")) ;

    ASSERT (Pattern) ;

    /* ---------------------------------------------------------------------- */
    /* print the sparse part of U */
    /* ---------------------------------------------------------------------- */

    k1 = 12 ;

    deg = Numeric->ulen ;
    if (deg > 0)
    {
	/* make last pivot row of U (singular matrices only) */
	for (j = 0 ; j < deg ; j++)
	{
	    Pattern [j] = Numeric->Upattern [j] ;
	}
    }

    PRINTF4 (("\n    row "ID":  length "ID".  End of Uchain.\n", INDEX (npiv-1),
	deg)) ;

    for (k = npiv-1 ; k >= 0 ; k--)
    {

	/* ------------------------------------------------------------------ */
	/* print row k of U */
	/* ------------------------------------------------------------------ */

	/* if prl is 3, print the first 10 entries of the first 10 columns */
	if (k1 > 0)
	{
	    prl = prl1 ;
	}

	up = Uip [k] ;
	ulen = Uilen [k] ;
	if (ulen < 0)
	{
	    return (FALSE) ;
	}
	newUchain = (up < 0) ;
	if (newUchain)
	{
	    up = -up ;
	    p = up + UNITS (Int, ulen) ;
	}
	else
	{
	    p = up ;
	}
	xp = (Entry *) (Numeric->Memory + p) ;
	if (deg > 0 && (p + (Int) UNITS (Entry, deg) > Numeric->size))
	{
	    return (FALSE) ;
	}
	for (j = 0 ; j < deg ; j++)
	{
	    col = Pattern [j] ;
	    PRINTF4 (("\tcol "ID" :", INDEX (col))) ;
	    if (prl >= 4) PRINT_ENTRY (*xp) ;
	    if (col <= k || col >= n_col)
	    {
		return (FALSE) ;
	    }
	    PRINTF4 (("\n")) ;
	    xp++ ;
	    /* truncate printout, but continue to check U */
	    if (prl == 4 && j == 9 && deg > 10)
	    {
		PRINTF (("\t...\n")) ;
		prl-- ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* make row k-1 of U in Pattern [0..deg-1] */
	/* ------------------------------------------------------------------ */

	if (k1-- > 0)
	{
	    prl = prl1 ;
	}
	else if (prl == 4)
	{
	    PRINTF (("    ...\n")) ;
	    prl-- ;
	}

	if (k > 0)
	{
	    PRINTF4 (("\n    row "ID":  ", INDEX (k-1))) ;
	}

	if (newUchain)
	{
	    /* next row is a new Uchain */
	    if (k > 0)
	    {
		deg = ulen ;
		PRINTF4 (("length "ID".  End of Uchain.\n", deg)) ;
		if (up + (Int) UNITS (Int, ulen) > Numeric->size)
		{
		    return (FALSE) ;
		}
		ip = (Int *) (Numeric->Memory + up) ;
		for (j = 0 ; j < deg ; j++)
		{
		    Pattern [j] = *ip++ ;
		}
	    }
	}
	else
	{
	    if (ulen > 0)
	    {
		PRINTF4 (("remove "ID" entries.  ", ulen)) ;
	    }
	    deg -= ulen ;
	    if (deg < 0)
	    {
		return (FALSE) ;
	    }
	    pos = Upos [k] ;
	    if (pos != EMPTY)
	    {
		/* add the pivot column */
		PRINTF4 (("add column "ID" at position "ID".  ",
		    INDEX (k), INDEX (pos))) ;
		if (pos < 0 || pos > deg)
		{
		    return (FALSE) ;
		}
		Pattern [deg++] = Pattern [pos] ;
		Pattern [pos] = k ;
	    }
	    PRINTF4 (("length "ID".\n", deg)) ;
	}
    }

    prl = prl1 ;
    PRINTF4 (("\n")) ;
    return (TRUE) ;
}

