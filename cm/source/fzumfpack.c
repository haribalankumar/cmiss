/*******************************************************************************
FILE : fzumfpack.c

LAST MODIFIED : 15 October 2002

DESCRIPTION : 

Set of stub routines for inclusion in cm when the Umfpack 4.0 library is
not available or included.
==============================================================================*/

/* Included files */

#include <stdio.h>

void umfpack_di_defaults(void)
{
	fprintf(stderr,"ERROR: link with Umfpack library: need umfpack_di_defaults\n");
}

int umfpack_di_symbolic(void)
{
	fprintf(stderr,"ERROR: link with Umfpack library: need umfpack_di_symbolic\n");
	return -1;
}

int umfpack_di_numeric(void)
{
	fprintf(stderr,"ERROR: link with Umfpack library: need umfpack_di_numeric\n");
	return -1;
}

int umfpack_di_solve(void)
{
	fprintf(stderr,"ERROR: link with Umfpack library: need umfpack_di_solve\n");
	return -1;
}

void umfpack_di_free_numeric(void)
{
	fprintf(stderr,"ERROR: link with Umfpack library: need umfpack_di_free_numeric\n");
}

void umfpack_di_free_symbolic(void)
{
	fprintf(stderr,"ERROR: link with Umfpack library: need umfpack_di_free_symbolic\n");
}
