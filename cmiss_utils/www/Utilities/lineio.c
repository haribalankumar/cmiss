/*
 *  Line IO.c
 *
 *    Line Input Module.
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: FILE, fopen(), flclose(), etc */

#include <errno.h>

#include <stdlib.h>
    /* For: malloc(), free(), NULL etc */

#include <string.h>
    /* For: strlen(), strcmp(), etc */

#include <time.h>
    /* For: time_t, time(), localtime(), strftime() etc */

#include "lineio.h"
    /* For: FILESTATE, getLine(), discardLine(), etc */



/* -- Private Method Prototypes -----------------------------------------------
*/

static char *defaultGetLine( FILESTATE *theState );
static void defaultDiscardLine( FILESTATE *theState );



/* -- Public Constructors and Destructors -------------------------------------
*/

FILESTATE *newFileState( char *fileName )
  {
  return newVirtualFileState( fileName, defaultGetLine, defaultDiscardLine );
  }


FILESTATE *newVirtualFileState( char *fileName, GETLINE_FUNCTION *getLine,
  DISCARDLINE_FUNCTION *discardLine )
  {
  FILESTATE *tempFileState;

  tempFileState = malloc( sizeof( FILESTATE ) );
  if( tempFileState == NULL )
    {
    fprintf( stderr, "memory allocation error in newVirtualFileState()\n" );
    return NULL;
    }

  tempFileState->maxLineLength = 512;
  tempFileState->currentLine = malloc( 512 );
  if( tempFileState->currentLine == NULL )
    {
    fprintf( stderr, "Memory allocation error in newVirtualFileState()\n" );
    free( tempFileState );
    return NULL;
    }

  tempFileState->theFile = fopen( fileName, "r" );
  if( tempFileState->theFile == NULL )
    {
    fprintf( stderr, "Failed to open file %s in newVirtualFileState(): %s\n", fileName, strerror(errno) );
    free( tempFileState->currentLine );
    free( tempFileState );
    return NULL;
    }

  tempFileState->currentFilePos = ftell( tempFileState->theFile );
      /* -- Unnessessary, see code below */

  tempFileState->readPending = True;

  tempFileState->getLine = getLine;
  tempFileState->discardLine = discardLine;

  return tempFileState;
  }


void disposeFileState( FILESTATE *theFileState )
  { 
  fclose( theFileState->theFile );
  free( theFileState->currentLine );
  free( theFileState );
  }



/* -- Private Methods ---------------------------------------------------------
*/

static char *defaultGetLine( FILESTATE *theState )
  {
  if( theState->readPending == True )
    {
    theState->currentFilePos = ftell( theState->theFile );
    fgets( theState->currentLine, theState->maxLineLength, theState->theFile );
    theState->readPending = False;
    }

  if( feof( theState->theFile ) )
    return NULL;

  return theState->currentLine;
  }


static void defaultDiscardLine( FILESTATE *theState )
  {
  theState->readPending = True;
  }



/* -- Public Methods ----------------------------------------------------------
*/ 

char *getLine( FILESTATE *theState )
  {
  return theState->getLine( theState );
  }


void discardLine( FILESTATE *theState )
  {
  theState->discardLine( theState );
  }


FILEPOS getFilePosition( FILESTATE *theState )
  {
  if( theState->readPending == True )
    (void) getLine( theState );

  return theState->currentFilePos;
  }


void setFilePosition( FILESTATE *theState, FILEPOS thePosition )
  {
  fseek( theState->theFile, thePosition, SEEK_SET );
  theState->currentFilePos = thePosition;
  theState->readPending = True;
  }
