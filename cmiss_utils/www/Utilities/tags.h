/*
 *  Tags.h
 *  
 *    A module for obtaining multiline tags from files
 */
#ifndef TAGS_H
#define TAGS_H


/* -- Required Includes -------------------------------------------------------
*/

#include "lineio.h"
    /* For: FILESTATE below */



/* -- Public Datatypes --------------------------------------------------------
*/

#if !defined( BOOLEAN )
#define BOOLEAN int
#define True 1
#define False 0
#endif


typedef struct
  {
  char *tag;
  char *body;
  }
TAG;


typedef struct TAGITEM
  {
  TAG *tag;

  struct TAGITEM *nextTagItem;
  }
TAGITEM;


typedef struct
  {
  int      numEntries; 
  TAGITEM *firstTagItem;
  TAGITEM *lastTagItem;
  }
TAGLIST;


typedef struct TAGFILE
  {
  TAG       *currentTag;
  FILEPOS    currentFilePos;
  FILESTATE *file;
  BOOLEAN    deallocateFileState; /* dispose FILESTATE, or leave for others? */
  BOOLEAN    readPending;
  }
TAGFILE;



/* -- Public Constructors, Destructors & Duplicators --------------------------
*/

TAG *newTag( void );
TAG *newTagFromContents( const char *name, const char *body );
void disposeTag( TAG *theTag );
TAG *duplicateTag( TAG *originalTag );

TAGITEM *newTagItem( TAG *tag );
void disopseTagItem( TAGITEM *tagItem );
void disposeTagItems( TAGITEM *tagItem );  /* recursively dispose items */

TAGLIST *newTagList( void );
void disposeTagList( TAGLIST *tagList );
TAGLIST *duplicateTagList( TAGLIST *tagList );

TAGFILE *newTagFile( char *fileName );
TAGFILE *newTagFileFromFileState( FILESTATE *fileState );
void disposeTagFile( TAGFILE *tagFileState );



/* -- Public Module Methods ---------------------------------------------------
*/

TAG *getTag( TAGFILE *tagFileState );
void discardTag( TAGFILE *tagFileState );

FILEPOS getTagFilePosition( TAGFILE *theFile );
void setTagFilePosition( TAGFILE *theFile, FILEPOS thePosition );

void printTag( TAG *theTag );
BOOLEAN tagCompare( const char *firstString, const char *secondString );

void addTagToList( TAGLIST *tagList, TAG *tag );


#endif
