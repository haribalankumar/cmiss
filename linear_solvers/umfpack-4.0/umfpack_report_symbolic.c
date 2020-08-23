/* ========================================================================== */
/* === UMFPACK_report_symbolic ============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Prints the Symbolic object. See umfpack_report_symbolic.h
    for details.

    Dynamic memory usage:  Allocates a size MAX (n_row,n_col)*sizeof(Int)
    workspace via a single call to UMF_malloc and then frees all of it via
    UMF_free on return.  The workspace is not allocated if an early error
    return occurs  before the workspace is needed.

*/

#include "umf_internal.h"
#include "umf_valid_symbolic.h"
#include "umf_report_perm.h"
#include "umf_malloc.h"
#include "umf_free.h"

GLOBAL Int UMFPACK_report_symbolic
(
    void *SymbolicHandle,
    const double Control [UMFPACK_CONTROL]
)
{
    Int n_row, n_col, nz, nchains, nfr, maxfrsize, maxnrows, maxncols, prl,
	k, chain, frontid, frontid1, frontid2, kk, *Chain_start, *W,
	*Chain_maxrows, *Chain_maxcols, *Front_npivcol, *Front_1strow,
	*Front_leftmostdesc, *Front_parent, done, status1, status2 ;
    SymbolicType *Symbolic ;

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

    PRINTF (("Symbolic object: ")) ;

    Symbolic = (SymbolicType *) SymbolicHandle ;
    if (!UMF_valid_symbolic (Symbolic))
    {
	PRINTF (("ERROR: invalid\n")) ;
	return (UMFPACK_ERROR_invalid_Symbolic_object) ;
    }

    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;

    nz = Symbolic->nz ;

    nchains = Symbolic->nchains ;
    nfr = Symbolic->nfr ;
    maxfrsize = Symbolic->maxfrsize ;
    maxnrows = Symbolic->maxnrows ;
    maxncols = Symbolic->maxncols ;

    Chain_start = Symbolic->Chain_start ;
    Chain_maxrows = Symbolic->Chain_maxrows ;
    Chain_maxcols = Symbolic->Chain_maxcols ;
    Front_npivcol = Symbolic->Front_npivcol ;
    Front_1strow = Symbolic->Front_1strow ;
    Front_leftmostdesc = Symbolic->Front_leftmostdesc ;
    Front_parent = Symbolic->Front_parent ;

    PRINTF4 (("\n    matrix to be factorized:\n")) ;
    PRINTF4 (("\tn_row: "ID" n_col: "ID"\n", n_row, n_col)) ;
    PRINTF4 (("\tnumber of entries: "ID"\n", nz)) ;
    PRINTF4 (("\tdense row contol parameter used: %g\n", Symbolic->drow)) ;
    PRINTF4 (("\tdense column contol parameter used: %g\n", Symbolic->dcol)) ;
    PRINTF4 (("\tblock size used for dense matrix kernels:  "ID"\n",
	Symbolic->nb)) ;

    PRINTF4 (("    variable-size part of Numeric object:\n")) ;
    PRINTF4 (("\tminimum initial size (Units): "ID"\n",
	Symbolic->num_mem_init_usage)) ;
    PRINTF4 (("\testimated peak size (Units):  %g\n",
	Symbolic->num_mem_usage_est)) ;
    PRINTF4 (("\testimated final size (Units): %g\n",
	Symbolic->num_mem_size_est)) ;
    if (Symbolic->num_mem_usage_est > Int_MAX / sizeof (Int))
    {
	PRINTF4 (("\tWarning: peak size (in bytes) exceeds maximum integer\n"));
    }
    PRINTF4 (("    peak memory usage for symbolic factorization (Units): %g\n",
	Symbolic->peak_sym_usage)) ;

    PRINTF4 (("    frontal matrices / supercolumns:\n")) ;
    PRINTF4 (("\tnumber of frontal chains: "ID"\n", nchains)) ;
    PRINTF4 (("\tnumber of frontal matrices: "ID"\n", nfr)) ;
    PRINTF4 (("\tlargest frontal matrix size (entries): "ID"\n", maxfrsize)) ;
    PRINTF4 (("\tlargest frontal matrix row dimension: "ID"\n", maxnrows)) ;
    PRINTF4 (("\tlargest frontal matrix column dimension: "ID"\n", maxncols)) ;

    k = 0 ;
    done = FALSE ;

    for (chain = 0 ; chain < nchains ; chain++)
    {
	frontid1 = Chain_start [chain] ;
	frontid2 = Chain_start [chain+1] - 1 ;
	PRINTF4 (("\n    Frontal chain: "ID".  Frontal matrices "ID" to "ID"\n",
	    INDEX (chain), INDEX (frontid1), INDEX (frontid2))) ;
	PRINTF4 (("\tLargest frontal matrix in Frontal chain: "ID"-by-"ID"\n",
	    Chain_maxrows [chain], Chain_maxcols [chain])) ;
	for (frontid = frontid1 ; frontid <= frontid2 ; frontid++)
	{
	    kk = Front_npivcol [frontid] ;
	    PRINTF4 (("\tFront: "ID"  pivot cols: "ID" (pivot columns "ID" to "
		ID")\n", INDEX (frontid), kk, INDEX (k), INDEX (k+kk-1))) ;
	    PRINTF4 (("\t    pivot row candidates: "ID" to "ID"\n",
		INDEX (Front_1strow [Front_leftmostdesc [frontid]]),
		INDEX (Front_1strow [frontid+1]-1))) ;
	    PRINTF4 (("\t    leftmost descendant: "ID"\n",
		INDEX (Front_leftmostdesc [frontid]))) ;
	    PRINTF4 (("\t    1st new candidate row : "ID"\n",
		INDEX (Front_1strow [frontid]))) ;
	    PRINTF4 (("\t    parent:")) ;
	    if (Front_parent [frontid] == EMPTY)
	    {
		PRINTF4 ((" (none)\n")) ;
	    }
	    else
	    {
		PRINTF4 ((" "ID"\n", INDEX (Front_parent [frontid]))) ;
	    }
	    done = (frontid == 20 && frontid < nfr-1 && prl == 4) ;
	    if (done)
	    {
		PRINTF4 (("\t...\n")) ;
		break ;
	    }
	    k += kk ;
	}
	if (Front_npivcol [nfr] != 0)
	{
	    PRINTF4 (("\tFront: "ID" placeholder for "ID" empty columns\n",
		INDEX (nfr), Front_npivcol [nfr])) ;
	}
	if (done)
	{
	    break ;
	}
    }

    W = (Int *) UMF_malloc (MAX (n_row, n_col), sizeof (Int)) ;
    if (!W)
    {
	PRINTF (("ERROR: out of memory to check Symbolic object\n\n")) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }

    PRINTF4 (("\nInitial column permutation, Qtree: ")) ;
    status1 = UMF_report_perm (n_col, Symbolic->Cperm_init, W, prl, 0) ;

    PRINTF4 (("\nInitial row permutation, Ptree: ")) ;
    status2 = UMF_report_perm (n_row, Symbolic->Rperm_init, W, prl, 0) ;

    (void) UMF_free ((void *) W) ;

    if (status1 != UMFPACK_OK || status2 != UMFPACK_OK)
    {
	return (UMFPACK_ERROR_invalid_Symbolic_object) ;
    }

    PRINTF4 (("    Symbolic object:  ")) ;
    PRINTF (("OK\n\n")) ;
    return (UMFPACK_OK) ;
}

