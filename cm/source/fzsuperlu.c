/*******************************************************************************
FILE : fzsuperlu.c

LAST MODIFIED : 15 October 2002

DESCRIPTION : 

Set of stub routines for inclusion in cm when the SuperLU library is
not available or included.
==============================================================================*/

/* Included files */

#include <stdio.h>

#if defined (unix) || defined (_AIX) || defined (WIN32)
/* name mappings */
#define Xerbla xerbla_
#endif /* defined (unix) */

void StatAlloc(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need StatAlloc\n");
}

void StatInit(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need StatInit\n");
}

void StatFree(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need StatFree\n");
}

void *intMalloc(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need intMalloc\n");
	return NULL;
}

void *superlu_malloc(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need superlu_malloc\n");
	return NULL;
}

void superlu_free(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need superlu_free\n");
}

void dCreate_Dense_Matrix(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need dCreate_Dense_Matrix\n");
}

void dCreate_CompCol_Matrix(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need dCreate_CompCol_Matrix\n");
}

void Destroy_CompCol_Matrix(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need Destroy_CompCol_Matrix\n");
}

void Destroy_CompCol_Permuted(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need Destroy_CompCol_Permuted\n");
}

void Destroy_CompCol_NCP(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need Destroy_CompCol_NCP\n");
}

void Destroy_SuperNode_Matrix(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need Destroy_SuperNode_Matrix\n");
}

void Destroy_SuperNode_SCP(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need Destroy_SuperNode_SCP\n");
}

void Destroy_SuperMatrix_Store(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need Destroy_SuperMatrix_Store\n");
}

void pdgstrf_init(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need pdgstrf_init\n");
}

void pdgstrf(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need pdgstrf\n");
}

void dgstrf(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need dgstrf\n");
}

void dgstrs(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need dgstrs\n");
}

void Xerbla_(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need xerbla\n");
}

void get_perm_c(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need get_perm_c\n");
}

void sp_preorder(void)
{
	fprintf(stderr,"ERROR: link with SuperLU library: need sp_preorder\n");
}
