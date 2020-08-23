/*
 *  Test Line IO.c
 *
 *    Test Line Input Module.
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fgets(), feof() etc */

#include <string.h>
    /* For: strlen(), strcmp(), etc */

#include "lineio.h"
    /* For: FILESTATE, newVirtualFileState() etc */



/* -- Private Method Prototypes -----------------------------------------------
*/

static BOOLEAN commentLine( char *line );
static char *testGetLine( FILESTATE *theState );
static void testDiscardLine( FILESTATE *theState );


/* -- Public Constructors and Destructors -------------------------------------
*/

FILESTATE *newTestFileState( char *fileName )
  {
  return newVirtualFileState( fileName, testGetLine, testDiscardLine );
  }



/* -- Private Methods ---------------------------------------------------------
*/

static BOOLEAN commentLine( char *line )
  {
  char *commentPrefix = "--";
  int   i;

  for( i = 0; i < 2; i++ )
    {
    if( line[i] == '\0' )
      return False;
    else
      if( line[i] != commentPrefix[i] )
         return False;
    }

  return True;
  }


static char *testGetLine( FILESTATE *theState )
  {
  BOOLEAN commentLineFound = True;

  if( theState->readPending == True )
    {
    while( commentLineFound == True )
      {
      fgets( theState->currentLine, theState->maxLineLength, theState->theFile);

      if( feof( theState->theFile ) )
        return NULL;
      
      if( commentLine( theState->currentLine ) == False )
        commentLineFound = False;
      }

    theState->readPending = False;
    }

  return theState->currentLine;
  }


static void testDiscardLine( FILESTATE *theState )
  {
  theState->readPending = True;
  }

