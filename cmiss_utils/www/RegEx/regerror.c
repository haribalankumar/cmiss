/*
 *  RegError.c
 * 
 *    Simple error reporting.
 *    Copyright Henry Spencer, ANSI-fication for DES by Edouard Poor.
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stderr */

#include <stdlib.h>
    /* For: exit() */
    


/* -- Module Public Methods ---------------------------------------------------
*/

void regerror( char *s )
  {
#ifdef ERRAVAIL
  error("regexp: %s", s);
#else
  fprintf(stderr, "regexp(3): %s", s);
#if EXIT_UPON_ERROR
  exit(1);
#endif
#endif
  /* NOTREACHED */
  }
