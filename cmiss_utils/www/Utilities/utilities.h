/*
 *  Utilities.h
 */

#ifndef UTILITIES_H
#define UTILITIES_H


/* -- Required Include Files --------------------------------------------------
*/

#include <stdio.h>
    /* For: FILE */



/* -- Public Types ------------------------------------------------------------
*/

#if !defined( BOOLEAN )
#define BOOLEAN int
#define True 1
#define False 0
#endif



/* -- Public Functions --------------------------------------------------------
*/ 

void forceLower( char *string );
BOOLEAN htmlString( char *string );
void outputString( FILE *output, char *string );


#endif
