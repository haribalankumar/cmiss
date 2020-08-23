#include <stdlib.h>
#include <stdio.h>

void *fmalloc_(long *n, long *nsize, int *init_flag)
{
  void *ptr;
  size_t nbyte;

  nbyte = (size_t) *n * (size_t) *nsize;
  /* fprintf(stderr,"falloc: allocating %ld bytes\n", nbyte); */

  if (*init_flag == 1) {
    ptr = calloc((size_t) *n, (size_t) *nsize);
  }
  else {
    ptr = malloc(nbyte);
  }

  return ptr;
}

void ffree_(void **ptr)
{
  free(*ptr);
}
