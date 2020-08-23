/*
 *  Print Error.c
 *
 *    Write out an error dialog
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: printf() */

#include <stdlib.h>
    /* For: exit() */

#include "print-error.h"
    /* For: Our public interface */



/* -- Module Datatypes --------------------------------------------------------
*/


/* -- Module Constants --------------------------------------------------------
*/ 


/* -- Public Modules --------------------------------------------------------
*/ 


void printErrorMessage( char *title, char *message )
  {
	fprintf( stdout,   "Content-type: text/html\n\n"
    "<HTML>\n<HEAD>\n<TITLE>%s</TITLE>\n</HEAD>\n"
    "<BODY BGCOLOR=\"#ffffff\">\n"
    "<H1>%s</H1>\n"
    "<P>%s<P></BODY></HTML>", title, title, message );
  fflush( stdout );
	exit(0);
  }


