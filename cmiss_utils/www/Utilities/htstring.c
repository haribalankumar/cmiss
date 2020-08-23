/*
 *  HTStrings.c
 *
 *    A module for handling
 *      o  x-url-encoded strings
 *      o  escaping strings to HTML-safe
 */


/* -- Required Includes -------------------------------------------------------
*/

#include <stdlib.h>
    /* For: malloc(), free(), NULL etc */

#include <stdio.h>
    /* For: printf() */

#include <string.h>
    /* For: strlen() etc */

#include "htstring.h"
    /* For: our public interface */



/* -- Macros and Defines ------------------------------------------------------
*/

#define EnsureSafetyTableInitialised() \
  if( urlEncodedTableInitialised == 0 ) \
    { \
    initEncodingSafetyTable( encodingSafetyTable, safeChars ); \
    urlEncodedTableInitialised = 1; \
    }

#define EnsureEscapeTableInitialised() \
  if( htmlEscapeTableInitialised == 0 ) \
    { \
    initEscapeTable( &htmlEscapeTable, htmlEscapeList ); \
    htmlEscapeTableInitialised = 1; \
    }



/* -- Private Datatypes -------------------------------------------------------
*/ 

#if !defined( BOOLEAN )
#define BOOLEAN int
#define True 1
#define False 0
#endif

/* Safety Table Stuff */

typedef enum
  {
  SafeEncoding = 1,
  EscapedEncoding = 2,
  PlusToSpaceEncoding = 3
  }
SAFETYENCODING;


/* Escape Table Stuff */

typedef enum
  {
  Safe = 1,
  Unsafe
  }
SAFE;

typedef struct
  {
  char  character;
  char *string;
  }
ESCAPELIST;

typedef struct
  {
  BOOLEAN     initalised;
  SAFE        safeTable[256];          /* UCHAR_MAX+1 instead of 256? */
  char       *replacementTable[256];
  int         lengthTable[256];
  }
ESCAPETABLE;



/* -- Private Module Variables ------------------------------------------------
*/

/* Safety Table Stuff */

static SAFETYENCODING encodingSafetyTable[256];
static int urlEncodedTableInitialised = 0;
static char *safeChars = "abcdefghijklmnopqrstuvwxyz"
                         "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                         "0123456789"
                         "-_@.~`'";

/* Escape Table Stuff */

static ESCAPETABLE htmlEscapeTable;
static int htmlEscapeTableInitialised = 0;
ESCAPELIST htmlEscapeList[] = { { '<',  "&lt;" },
                                { '>',  "&gt;" },
                                { '&',  "&amp;" },
                                { '\"', "&quot;" },
                                { 0,    NULL } }; 



/* -- Private Methods ---------------------------------------------------------
*/

void initEncodingSafetyTable( SAFETYENCODING theTable[256], char *safeChars )
  {
  int i;

  for( i = 0; i < 256; i++ )
    {
    theTable[i] = EscapedEncoding;
    }

  for( i = 0; i < strlen( safeChars ); i++ )
    {
    theTable[ (unsigned int) safeChars[i] ] = SafeEncoding;
    }

  theTable[ (unsigned int) ' ' ] = PlusToSpaceEncoding;
  }


void initEscapeTable( ESCAPETABLE *theTable, ESCAPELIST escapeList[] )
  {
  int i;
  unsigned int charValue;

  for( i = 0; i < 256; i++ )
    {
    theTable->safeTable[i] = Safe;
    theTable->replacementTable[i] = NULL;
    theTable->lengthTable[i] = 1;
    }

  i = 0;
  while( escapeList[i].string != NULL )
    {
    charValue = (unsigned int) escapeList[i].character;

    theTable->safeTable[ charValue ] = Unsafe;
    theTable->replacementTable[ charValue ] = escapeList[i].string;
    theTable->lengthTable[ charValue ] = strlen( escapeList[i].string );
    i++;
    }
  }



/* -- Public Methods for urlEncoded[...] --------------------------------------
*/

size_t urlEncodedStrlen( const char *theString )
  {
  int         length;
  const char *ptr;

  EnsureSafetyTableInitialised();

  length = 0;
  ptr = theString;
  while( *ptr != '\0' )
    {
    switch( encodingSafetyTable[ (unsigned int) *ptr ] )
      {
      case SafeEncoding:
        length += 1;
        break;

      case EscapedEncoding:
        length += 3;
        break;

      case PlusToSpaceEncoding:
        length += 1;
        break;

      default:
        fprintf( stderr, "Error in urlEncodedStrlen() - unknown value\n" );
      }
   
    ptr++;
    }

  return length;
  }


char *urlEncodedStrcpy( char *destinationString, const char *sourceString )
  {
  char *destPtr;
  const char *sourcePtr;
  char hexTable[] = "0123456789ABCDEF"; 

  EnsureSafetyTableInitialised();

  destPtr = destinationString;
  sourcePtr = sourceString;

  while( *sourcePtr != '\0' )
    {
    switch( encodingSafetyTable[ (unsigned int) *sourcePtr ] )
      {
      case SafeEncoding:
        *destPtr++ = *sourcePtr++;
        break;
      
      case EscapedEncoding: 
        *destPtr++ = '%';
        *destPtr++ = hexTable[ (unsigned int) (((*sourcePtr) / 16) % 16) ];
        *destPtr++ = hexTable[ (unsigned int) ((*sourcePtr++) % 16) ];
          /* -- the 16's above are a bit 8-bit char-o-centric... change? */
        break;

      case PlusToSpaceEncoding:
        *destPtr++ = '+';
        sourcePtr++;
        break;

      default:
        fprintf( stderr, "Error in urlEncodedStrcpy() - unknown value\n" );
      }
    }
  *destPtr = '\0';

  return destinationString;
  }


char *urlEncodedStrcat( char *destinationString, const char *sourceString )
  {
  char *endPtr;

  endPtr = destinationString + strlen( destinationString );
  
  urlEncodedStrcpy( endPtr, sourceString );

  return destinationString;
  }



/* -- Public Methods for htmlEscaped[...] -------------------------------------
*/

size_t htmlEscapedStrlen( const char *sourceString )
  {
  const char *sourcePtr;
  size_t      length;

  EnsureEscapeTableInitialised();

  sourcePtr = sourceString;
  length = 0;
  while( *sourcePtr != '\0' )
    {
    length += htmlEscapeTable.lengthTable[ (unsigned int) *sourcePtr ];
    sourcePtr++;
    }

  return length;
  }


char *htmlEscapedStrcpy( char *destinationString, const char *sourceString )
  {
  char       *destPtr;
  const char *sourcePtr;
  char       *replacementPtr;

  EnsureEscapeTableInitialised();

  destPtr = destinationString;
  sourcePtr = sourceString;

  while( *sourcePtr != '\0' )
    {
    if( htmlEscapeTable.safeTable[ *sourcePtr ] == Safe )
      {
      *destPtr++ = *sourcePtr++;
      }
    else
      {
      replacementPtr = htmlEscapeTable.replacementTable[ *sourcePtr ];

      while( *replacementPtr != '\0' )
        {
        *destPtr++ = *replacementPtr++;
        }

      sourcePtr++;
      }
    }
  *destPtr = '\0';

	return (destinationString);
  }


char *htmlEscapedStrcat( char *destinationString, const char *sourceString )
  {
  char *endPtr;

  endPtr = destinationString + strlen( destinationString );

  htmlEscapedStrcpy( endPtr, sourceString );

  return destinationString;
  }

