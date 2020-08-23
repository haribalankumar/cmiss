/* ========================================================================== */
/* === UMF_report_vector ==================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

#include "umf_internal.h"

GLOBAL Int UMF_report_vector
(
    Int n,
    const double Xx [ ],
    const double Xz [ ],
    Int prl,
    Int user
)
{
    Int n2, i ;
    Entry *X, xi ;

    /* if Xz is null, then X is in "merged" format (compatible with Entry, */
    /* and ANSI C99 double _Complex type). */
    X = (Entry *) Xx ;

    ASSERT (prl >= 3) ;

    if (user || prl >= 4)
    {
	PRINTF (("dense vector, n = "ID". ", n)) ;
    }

    if (user)
    {
	if (!Xx)
	{
	    PRINTF (("ERROR: vector not present\n\n")) ;
	    return (UMFPACK_ERROR_argument_missing) ;
	}
	if (n < 0)
	{
	    PRINTF (("ERROR: length of vector is < 0\n\n")) ;
	    return (UMFPACK_ERROR_n_nonpositive) ;
	}
    }

    if (user || prl >= 4)
    {
	PRINTF4 (("\n")) ;
    }

    if (prl == 4)
    {
	n2 = MIN (10, n) ;
	for (i = 0 ; i < n2 ; i++)
	{
	    PRINTF (("    "ID" :", INDEX (i))) ;
	    if (Xz)
	    {
	    	ASSIGN (xi, Xx [i], Xz [i]) ;
	    }
	    else
	    {
	    	xi = X [i] ;
	    }
	    PRINT_ENTRY (xi) ;
	    PRINTF (("\n")) ;
	}
	if (n2 < n)
	{
	    PRINTF (("    ...\n")) ;
	    PRINTF (("    "ID" :", INDEX (n-1))) ;
	    if (Xz)
	    {
	    	ASSIGN (xi, Xx [n-1], Xz [n-1]) ;
	    }
	    else
	    {
	    	xi = X [n-1] ;
	    }
	    PRINT_ENTRY (xi) ;
	    PRINTF (("\n")) ;
	}
    }
    else if (prl > 4)
    {
	/* print level 4 or more */
	for (i = 0 ; i < n ; i++)
	{
	    PRINTF (("    "ID" :", INDEX (i))) ;
	    if (Xz)
	    {
	    	ASSIGN (xi, Xx [i], Xz [i]) ;
	    }
	    else
	    {
	    	xi = X [i] ;
	    }
	    PRINT_ENTRY (xi) ;
	    PRINTF (("\n")) ;
	}
    }

    PRINTF4 (("    dense vector ")) ;
    if (user || prl >= 4)
    {
	PRINTF (("OK\n\n")) ;
    }
    return (UMFPACK_OK) ;
}

