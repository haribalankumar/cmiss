/*
 *  Line IO.h
 *
 *    Public interface to the Line IO module
 */
#ifndef LINEIO_H
#define LINEIO_H


/* -- Required Include Files --------------------------------------------------
*/

#include <stdio.h>
    /* For: FILE below */



/* -- Public Datatypes --------------------------------------------------------
*/

#if !defined( BOOLEAN )
#define BOOLEAN int
#define True 1
#define False 0
#endif

struct FILESTATE;

typedef char *GETLINE_FUNCTION( struct FILESTATE * );
typedef void DISCARDLINE_FUNCTION( struct FILESTATE * );

typedef long FILEPOS;

typedef struct FILESTATE
  {
  /* Private Data */
  FILE    *theFile;
  FILEPOS  currentFilePos;
  char    *currentLine;
  int      maxLineLength;
  BOOLEAN  readPending;

  /* Private Virtual Calls */
  GETLINE_FUNCTION     *getLine;
  DISCARDLINE_FUNCTION *discardLine;
  }
FILESTATE;



/* -- Constructors and Destructors --------------------------------------------
*/

FILESTATE *newFileState( char *fileName );
FILESTATE *newVirtualFileState( char *fileName, GETLINE_FUNCTION *getLine,
  DISCARDLINE_FUNCTION *discardLine );
void disposeFileState( FILESTATE *theFileState );



/* -- Module Method Prototypes ------------------------------------------------
*/

char *getLine( FILESTATE *theState );
void discardLine( FILESTATE *theState );
FILEPOS getFilePosition( FILESTATE *theState );
void setFilePosition( FILESTATE *theState, FILEPOS filePos );


#endif
