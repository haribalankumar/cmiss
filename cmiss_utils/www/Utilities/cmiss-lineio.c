/*
 *  Line IO.c
 *
 *    Line Input Module.
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
static char *cmissGetLine( FILESTATE *theState );
static void cmissDiscardLine( FILESTATE *theState );


/* -- Public Constructors and Destructors -------------------------------------
*/

FILESTATE *newCmissFileState( char *fileName )
  {
  return newVirtualFileState( fileName, cmissGetLine, cmissDiscardLine );
  }



/* -- Private Methods ---------------------------------------------------------
*/

static BOOLEAN commentLine( char *line )
  {
  char *commentPrefixOne = "C###  ";
  char *commentPrefixTwo = "c#### ";
  int   i;

  for( i = 0; i < 6; i++ )
    {
    if( line[i] == '\0' )
      return False;
    else
      if( line[i] != commentPrefixOne[i] && line[i] != commentPrefixTwo[i] )
         return False;
    }

  return True;
  }


static char *cmissGetLine( FILESTATE *theState )
  {
  BOOLEAN commentLineFound = False;

  if( theState->readPending == True )
    {
    while( commentLineFound == False )
      {
      fgets( theState->currentLine, theState->maxLineLength, theState->theFile);

      if( feof( theState->theFile ) )
        return NULL;
      
      if( commentLine( theState->currentLine ) == True )
        commentLineFound = True;
      }
    theState->readPending = False;
    }

  return theState->currentLine + 6;
  }


static void cmissDiscardLine( FILESTATE *theState )
  {
  theState->readPending = True;
  }

