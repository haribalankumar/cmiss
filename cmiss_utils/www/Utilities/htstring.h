/*
 *  HTStrings.h
 *
 *    Public interface to the HTStrings module
 */
#ifndef HTSTRING_H
#define HTSTRING_H


/* -- Required Includes -------------------------------------------------------
*/

#include <stddef.h>
    /* For: size_t used below */



/* -- Public Methods ----------------------------------------------------------
*/

size_t urlEncodedStrlen( const char *theString );
char *urlEncodedStrcpy( char *destinationString, const char *sourceString );
char *urlEncodedStrcat( char *destinationString, const char *sourceString );

size_t htmlEscapedStrlen( const char *sourceString );
char *htmlEscapedStrcpy( char *destinationString, const char *sourceString );
char *htmlEscapedStrcat( char *destinationString, const char *sourceString );


#endif
