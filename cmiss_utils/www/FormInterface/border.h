/*
 *  Border.h
 *
 *    Public Interface to the Border routines
 */
#ifndef BORDER_H
#define BORDER_H

#include <stdio.h>
    /* For: FILE* */

#include "user-managment.h"
    /* For: USER* */

/* -- Public Methods ----------------------------------------------------------
*/

void writeError(char *error_str);
void borderWriteError( FILE *output, USER *user, char *error_str );
int borderGenerateHeader( FILE *output, USER *user, BOOLEAN refresh,
                          char *address, char *other );
int borderGenerateFooter( FILE *output );

int paramBorderGenerateHeader( FILE *output, USER *user, BOOLEAN refresh,
                               char *address, char *other );
void paramBorderWriteError( FILE *output, USER *user, char *error_str );
#endif
