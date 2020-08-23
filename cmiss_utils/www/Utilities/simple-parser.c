/*
 *  Simple Parser.c
 *
 *    A very simle parsing "class".
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: FILE, fopen(), flclose(), etc */

#include <stdlib.h>
    /* For: malloc(), free(), NULL etc */

#include <string.h>
    /* For: strlen(), strcmp(), etc */

#include "simple-parser.h"
    /* For: Our public interface */



/* -- Private Method Prototypes -----------------------------------------------
*/

BOOLEAN listContainsChar( char *theList, char theChar );



/* -- Public Constructors/Destructors -----------------------------------------
*/

PARSESTATE *newParseState( char *theStream )
  {
  PARSESTATE *tempParseState;
  int         streamLength;

  tempParseState = malloc( sizeof( PARSESTATE ) );
  
  if( tempParseState == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newParseState() [1]\n" );
    return NULL;
    }

  streamLength = strlen( theStream );
  tempParseState->stream = malloc( streamLength + 1 );
  
  if( tempParseState->stream == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newParseState() [2]\n" );
    free( tempParseState );
    return NULL;
    }

  strcpy( tempParseState->stream, theStream );

  tempParseState->currentPtr = tempParseState->stream;
  tempParseState->end = tempParseState->stream + streamLength;

  return tempParseState;
  }


void disposeParseState( PARSESTATE *theParseState )
  {
  free( theParseState->stream );
  free( theParseState );
  }



/* -- Private Methods ---------------------------------------------------------
*/

BOOLEAN listContainsChar( char *theList, char theChar )
  {
  char *theListPtr = theList;

  while( *theListPtr != '\0' )
    if( theChar == *theListPtr )
      return True;
    else
      theListPtr++;

  /* We couldn't find it in the list, but a '\0' is implicitally in the list */
  if( theChar == '\0' )
    return True;
  else
    return False;
  }



/* -- Public Methods ----------------------------------------------------------
*/

BOOLEAN prefix( PARSESTATE *theState, char *thePrefix )
  {
  if( strncmp( theState->currentPtr, thePrefix, strlen(thePrefix) ) == 0 )
    return True;
  else
    return False;
  }


BOOLEAN empty( PARSESTATE *theState )
  {
  if( theState->currentPtr == theState->end )
    return True;
  else
    return False;
  }


void shiftChars( PARSESTATE *theState, int shiftAmount )
  {
  if( shiftAmount > theState->end - theState->currentPtr )
    theState->currentPtr = theState->end;
  else
    theState->currentPtr += shiftAmount;
  }


void skip( PARSESTATE *theState, char *skipList )
  {
  while( *theState->currentPtr != '\0' ) 
    {
    if( listContainsChar( skipList, *theState->currentPtr ) == False )
      return;
    else
      theState->currentPtr++;
    }
  }


void skipWord( PARSESTATE *theState, char *theWord )
  {
  int wordLength;

  wordLength = strlen( theWord );
  if( strncmp( theState->currentPtr, theWord, wordLength ) == 0 )
    theState->currentPtr += wordLength;
  }


char *extract( PARSESTATE *theState, char *terminationList )
  {
  char *returnString;
  int   stringLength = 0;
  char *searchPtr = theState->currentPtr;

  while( listContainsChar( terminationList, *searchPtr ) == False )
    {
    stringLength++;
    searchPtr++;
    }

  returnString = malloc( sizeof( char ) * (stringLength + 1) );

  strncpy( returnString, theState->currentPtr, stringLength );
  returnString[stringLength] = '\0';

  if( *searchPtr == '\0' )
    theState->currentPtr = searchPtr;
  else
    theState->currentPtr += stringLength; /* dont skip terminator -- ??? */

  return returnString;
  }


void skipUntil( PARSESTATE *theState, char *terminationList )
  {
  char *searchPtr = theState->currentPtr;

  while( listContainsChar( terminationList, *searchPtr ) == False )
    searchPtr++;

  theState->currentPtr = searchPtr;
  }

