/*
 *  Tags.c
 *  
 *    A module for obtaining multiline tags from files
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: puts() */

#include <errno.h>

#include <stdlib.h>
    /* For: NULL */

#include <string.h>
    /* For: strchr(), strdup(), strerror() */

#include <ctype.h>
    /* For: tolower() */

#include "simple-parser.h"
    /* For: PARSESTATE, newParseState(), prefix(), extract() etc */

#include "lineio.h"
    /* For: FILESTATE, getLine(), discardLine() etc */

#include "tags.h"
    /* For: our public interface */



/* -- Private Datatypes -------------------------------------------------------
*/

typedef enum
  {
  LegalTagChar,
  IllegalTagChar
  }
LEGALTAG;


typedef enum
  {
  TagWhitespace,
  TagNonWhitespace
  }
TAGWHITESPACE;


typedef struct BODYITEM
  {
  char            *text;
  int              textLength;
  struct BODYITEM *next;
  }
BODYITEM;


typedef struct 
  {
  size_t    totalLength;
  BODYITEM *first;
  BODYITEM *last;
  }
BODYLIST;



/* -- Private Constructors & Destructors --------------------------------------
*/ 

BODYITEM *newBodyItem( char *text );
void disposeBodyItem( BODYITEM *theBodyItem );
BODYLIST *newBodyList( void );
void disposeBodyListItems( BODYITEM *theItem );
void disposeBodyList( BODYLIST *theBodyList );
char *newStringFromBodyList( BODYLIST *bodyList );



/* -- Private Module Method Prototypes ----------------------------------------
*/

void initLegalTagTable( LEGALTAG table[256] );
void initWhitespaceTable( TAGWHITESPACE table[256] );
void skipInputForNextTag( FILESTATE *theState );
BOOLEAN getNextTag( FILESTATE *fileState, TAG *tag );
void addToBodyList( BODYLIST *bodyList, char *body );



/* -- Private Module Data -----------------------------------------------------
*/

static int legalTagTableInitialized = 0;
LEGALTAG legalTagTable[256];
char legalTagCharacters[] = "abcdefghijklmnopqrstuvwxyz"
                            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                            "1234567890"
                            "_-";

static int tagWhitespaceTableInitialized = 0;
TAGWHITESPACE tagWhitespaceTable[256]; 
char tagWhitespaceCharacters[] = "-_ \t\n\r";  /* Added "\n\r" Aug 19 1996 */



/* -- Public Constructors, Destructors & Duplicators --------------------------
*/

TAG *newTag( void )
  {
  TAG *tempTag;

  tempTag = malloc( sizeof( TAG ) );
  if( tempTag == NULL )
    {
    fprintf( stderr, "Memory allocation error in newTag()\n" );
    return NULL;
    }

  tempTag->tag = NULL;
  tempTag->body = NULL;
   
  return tempTag;
  }

#define FUNCTION "newTagFromContents"
TAG *newTagFromContents( const char *name, const char *body )
  {
  TAG *tempTag;

  tempTag = newTag();
  if( tempTag == NULL )
    {
    fprintf( stderr, "newTag() returned NULL in " FUNCTION "\n" );
    return NULL;
    }

  tempTag->tag = strdup( name );
  if( tempTag->tag == NULL )
    {
    fprintf( stderr, "Failed to duplicate name in " FUNCTION ": %s\n",
	     strerror(errno) );
    disposeTag( tempTag );
    return NULL;
    }

  tempTag->body = strdup( body );
  if( tempTag->body == NULL )
    {
    fprintf( stderr, "Failed to duplicate body in " FUNCTION ": %s\n",
	     strerror(errno) );
    disposeTag( tempTag );
    return NULL;
    }

  return tempTag;
  }
#undef FUNCTION

void disposeTag( TAG *theTag )
  {
  if( theTag->tag != NULL ) free( theTag->tag );
  if( theTag->body != NULL ) free( theTag->body );
  
  free( theTag );
  }


TAG *duplicateTag( TAG *originalTag )
  {
  TAG *copyTag;

  copyTag = newTag();
  if( copyTag == NULL )
    {
    fprintf( stderr, "newTag() returned NULL in duplicateTag()\n" );
    return NULL;
    }

  copyTag->tag = malloc( strlen( originalTag->tag ) + 1 );
  if( copyTag->tag == NULL )
    {
    fprintf( stderr, "Memory allocation failure in duplicateTag() [1]\n" );
    disposeTag( copyTag );
    return NULL;
    }
  
  strcpy( copyTag->tag, originalTag->tag );

  copyTag->body = malloc( strlen( originalTag->body ) + 1 );
  if( copyTag->body == NULL )
    {
    fprintf( stderr, "Memory allocation failure in duplicateTag() [2]\n" );
    disposeTag( copyTag );
    return NULL;
    }

  strcpy( copyTag->body, originalTag->body );

  return copyTag;
  }


TAGITEM *newTagItem( TAG *tag )
  {
  TAGITEM *tempTagItem;

  tempTagItem = malloc( sizeof( TAGITEM ) );
  if( tempTagItem == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newTagItem()\n" );
    return NULL;
    }

  tempTagItem->tag = tag;
  tempTagItem->nextTagItem = NULL;

  return tempTagItem;
  }


void disposeTagItem( TAGITEM *tagItem )
  {
  disposeTag( tagItem->tag );
  free( tagItem );
  }


void disposeTagItems( TAGITEM *tagItem )
  {
  if( tagItem != NULL )
    {
    disposeTagItems( tagItem->nextTagItem );
    disposeTagItem( tagItem );
    }
  }


TAGLIST *newTagList( void )
  {
  TAGLIST *tempTagList;

  tempTagList = malloc( sizeof( TAGLIST ) );
  if( tempTagList == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newTagList()\n" );
    return NULL;
    }

  tempTagList->numEntries = 0;
  tempTagList->firstTagItem = NULL;
  tempTagList->lastTagItem = NULL;

  return tempTagList;
  }


TAGLIST *duplicateTagList( TAGLIST *original )
  {
  TAGLIST *copy;
  TAG     *tagCopy;
  TAGITEM *tagItem;
  TAGITEM *tagItemCopy;

  copy = newTagList();
  if( copy == NULL )
    {
    fprintf( stderr, "newTagList() returned NULL in duplicateTagList()\n" );
    return NULL;
    }

  tagItem = original->firstTagItem;
  while( tagItem != NULL )
    {
    tagCopy = duplicateTag( tagItem->tag );
    if( tagCopy == NULL )
      {
      fprintf( stderr, "duplicateTag() returned NULL in "
        "duplicateTagList()\n" );
      disposeTagList( copy );
      return NULL;
      }

    tagItemCopy = newTagItem( tagCopy );
    if( tagItemCopy == NULL )
      {
      fprintf( stderr, "newTagItem() returned NULL in duplicateTagList()\n" );
      disposeTagList( copy );
      return NULL;
      }

    if( copy->lastTagItem != NULL )
      {
      copy->lastTagItem->nextTagItem = tagItemCopy;
      copy->lastTagItem = tagItemCopy;
      }
    else
      {
      copy->firstTagItem = copy->lastTagItem = tagItemCopy;
      }
    
    copy->numEntries++;

    tagItem = tagItem->nextTagItem;
    }

  return copy;
  }


void disposeTagList( TAGLIST *tagList )
  {
  disposeTagItems( tagList->firstTagItem );
  free( tagList );
  }


TAGFILE *newTagFileState( char *fileName )
  {
  TAGFILE *tempTagFileState;

  tempTagFileState = malloc( sizeof( TAGFILE ) );
  if( tempTagFileState == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newTagFileState()\n" );
    return NULL;
    }

  tempTagFileState->file = newFileState( fileName );
  if( tempTagFileState->file == NULL )
    {
    fprintf( stderr, "newFileState() return NULL in newTagFileState()" );
    free( tempTagFileState );
    return NULL;
    }

  tempTagFileState->currentFilePos = getFilePosition( tempTagFileState->file );
      /* -- Unnessessary, is overwritten later */

  tempTagFileState->readPending = True;
  tempTagFileState->currentTag = NULL;

  return tempTagFileState;
  }


TAGFILE *newTagFileStateFromFileState( FILESTATE *fileState )
  {
  TAGFILE *tempTagFileState;

  tempTagFileState = malloc( sizeof( TAGFILE ) );
  if( tempTagFileState == NULL )
    { 
    fprintf( stderr, "Memory allocation failure in newTagFileState()\n" );
    return NULL;
    }
 
  tempTagFileState->file = fileState;
  tempTagFileState->readPending = True;
  tempTagFileState->currentTag = NULL;

  return tempTagFileState;
  }


void disposeTagFileState( TAGFILE *tagFileState )
  {
  disposeFileState( tagFileState->file );
  if( tagFileState->currentTag != NULL )
    disposeTag( tagFileState->currentTag );

  free( tagFileState );
  }



/* -- Private Constructors & Destructors --------------------------------------
*/  

BODYITEM *newBodyItem( char *text ) 
  {
  BODYITEM *tempBodyItem;
  int       length;

  tempBodyItem = malloc( sizeof( BODYITEM ) );
  if( tempBodyItem == NULL )
    {
    fprintf( stderr, "Memory allocation error in newBodyItem\n" );
    return NULL;
    }

  length = strlen( text );

  tempBodyItem->text = malloc( length + 1 );
  if( tempBodyItem->text == NULL )
    {
    fprintf( stderr, "Memory allocation error in newBodyItem\n" );
    free( tempBodyItem );
    return NULL;
    }

  strcpy( tempBodyItem->text, text );
  tempBodyItem->textLength = length;
  tempBodyItem->next = NULL;

  return tempBodyItem;
  }


void disposeBodyItem( BODYITEM *theBodyItem )
  {
  if( theBodyItem->text != NULL )
    free( theBodyItem->text );

  free( theBodyItem );
  }


BODYLIST *newBodyList( void )
  {
  BODYLIST *tempBodyList;

  tempBodyList = malloc( sizeof( BODYLIST ) );
  if( tempBodyList == NULL )
    {
    fprintf( stderr, "Memory allocation error in newBodyList\n" );
    return NULL;
    }
 
  tempBodyList->first = NULL;
  tempBodyList->last = NULL;
  tempBodyList->totalLength = 0;

  return tempBodyList;
  }


void disposeBodyListItems( BODYITEM *theItem )
  {
  if( theItem != NULL )
    {
    disposeBodyListItems( theItem->next );
    disposeBodyItem( theItem );
    }
  }


void disposeBodyList( BODYLIST *theBodyList )
  {
  disposeBodyListItems( theBodyList->first );
  free( theBodyList );
  }


char *newStringFromBodyList( BODYLIST *bodyList )
  {
  char     *tempString;
  BODYITEM *bodyPtr;
  char     *stringPtr;

  tempString = malloc( bodyList->totalLength + 1 );
  if( tempString == NULL )
    {
    fprintf( stderr, "Memory allocation error in newStringFromBodyList()\n" );
    return NULL;
    }
  
  bodyPtr = bodyList->first;
  stringPtr = tempString;

  while( bodyPtr != NULL )
    {
    memcpy( stringPtr, bodyPtr->text, bodyPtr->textLength );
    stringPtr += bodyPtr->textLength;

    bodyPtr = bodyPtr->next;
    }

  *stringPtr = '\0';

  return tempString;
  }



/* -- Private Module Methods --------------------------------------------------
*/

/*
 *  Question: What is the type of a char constant? Is it signed or unsigned
 *  int? Can the char constant ever evaluate to a negative number?
 */
void initLegalTagTable( LEGALTAG table[256] )
  {
  unsigned int i;

  for( i = 0; i < 256; i++ )
    table[ i ] = IllegalTagChar;

  for( i = 0; i < strlen( legalTagCharacters ); i++ )
    table[ (unsigned int) legalTagCharacters[i] ] = LegalTagChar;
  }


void initWhitespaceTable( TAGWHITESPACE table[256] )
  {
  unsigned int i;

  for( i = 0; i < 256; i++ )
    table[ i ] = TagNonWhitespace;

  for( i = 0; i < strlen( tagWhitespaceCharacters ); i++ )
    table[ (unsigned int) tagWhitespaceCharacters[i] ] = TagWhitespace;
  }


void skipInputForNextTag( FILESTATE *theState )
  {
  char    *colonPos;
  BOOLEAN  tagLine;
  char    *theLine;
  char    *ptr;

  if( legalTagTableInitialized == 0 )
    {
    initLegalTagTable( legalTagTable );
    legalTagTableInitialized = 1;
    }

  theLine = getLine( theState );

  while( theLine != NULL )
    {
    colonPos = strchr( theLine, ':' );
    if( colonPos != NULL )
      {
      tagLine = True;  /* We then try to disprove this */

      for( ptr = theLine; ptr < colonPos; ptr++ )
        if( legalTagTable[ *ptr ] == IllegalTagChar )
          {
          tagLine = False;
          break;
          }

      if( tagLine == True )
        return;
      }

    discardLine( theState );

    theLine = getLine( theState );
    }
  }


/*
 *  This routine is in serious need of help.
 * 
 *  Well I rewrote this routine. It got worse. Sorry.
 */
BOOLEAN getNextTag( FILESTATE *fileState, TAG *tag )
  {
  char       *line;
  char       *ptr;
  int         length;
  char       *colonPos;
  char       *body;
  BODYLIST   *bodyList;
  BOOLEAN     firstLinePending;
  int         count;
  int         hangingIndent;

  line = getLine( fileState );
  if( line == NULL )
    return False;

  firstLinePending = True;
  hangingIndent = 0;

  /* == Search for "name: [content]" line == */

  bodyList = newBodyList();

  colonPos = strchr( line, ':' );
  if( colonPos == NULL )
    return False;

  length = colonPos - line;
  tag->tag = malloc( length + 1 );
  strncpy( tag->tag, line, length );
  tag->tag[length] = '\0';

  ptr = colonPos + 1;
  while( *ptr == ' ' || *ptr == '\t' )
    ptr++;

  if( *ptr != '\n' )
    {
    addToBodyList( bodyList, ptr );

    firstLinePending = False;
    hangingIndent = 999; /* Nasty large number */
    }

  discardLine( fileState );
    /* -- Discard the line */

  line = getLine( fileState );
    /* -- get next line */

  while( line != NULL )
    {
    if( line[0] == ' ' || line[0] == '\t' || line[0] == '\n' ) 
      {
      if( firstLinePending == True )
        {
        ptr = line;
        while( *ptr == ' ' || *ptr == '\t' )
          ptr++;

        if( *ptr != '\n' )
          {
          firstLinePending = False;
          hangingIndent = ptr - line;
     
          addToBodyList( bodyList, ptr );
          }
        }
      else
        {
        ptr = line;
        count = hangingIndent;
        while( ( count-- != 0 ) && ( *ptr == ' ' || *ptr == '\t' ) )
          ptr++;

        addToBodyList( bodyList, ptr );
        }

      discardLine( fileState );
      line = getLine( fileState );
      }
    else
      {
      break;
      }
    }

  body = newStringFromBodyList( bodyList );
  disposeBodyList( bodyList );

  tag->body = body;
  return True;
  }


void addToBodyList( BODYLIST *bodyList, char *body )
  {
  bodyList->totalLength += strlen( body );
  
  if( bodyList->first == NULL )
    {
    bodyList->first = newBodyItem( body );
    bodyList->first->next = NULL;
    bodyList->last = bodyList->first;
    }
  else
    {
    bodyList->last->next = newBodyItem( body );
    bodyList->last->next->next = NULL;
    bodyList->last = bodyList->last->next;
    }
  }



/* -- Public Module Methods ---------------------------------------------------
*/

TAG *getTag( TAGFILE *tagFile )
  {
  if( tagFile->readPending == True )
    {
    if( tagFile->currentTag != NULL )
      disposeTag( tagFile->currentTag );

    tagFile->currentTag = newTag();

    skipInputForNextTag( tagFile->file );

    tagFile->currentFilePos = getFilePosition( tagFile->file );

    if( getNextTag( tagFile->file, tagFile->currentTag ) == False )
      {
      /* 
       *  getNextTag() failed. (This is ugly - getNextTag() should be
       *  rewritten not to have these multiple return symantics.)
       */
      disposeTag( tagFile->currentTag );
      tagFile->currentTag = NULL;
      }

    tagFile->readPending = False;
    }

  return tagFile->currentTag;
  }


void discardTag( TAGFILE *tagFile )
  {
  tagFile->readPending = True;
  }


TAGFILE *newTagFile( char *fileName )
  {
  FILESTATE *fileState;
  TAGFILE   *tempTagFile;

  fileState = newFileState( fileName );
  if( fileState == NULL )
    {
    fprintf( stderr, "newFileState() returned NULL in newTagFile()\n" );
    return NULL;
    }

  tempTagFile = newTagFileFromFileState( fileState );
  if( tempTagFile == NULL )
    {
    fprintf( stderr, "newTagFileFromFileState() return NULL in "
      "newTagFile()\n" );
    disposeFileState( fileState );
    return NULL;
    }

  tempTagFile->deallocateFileState = True;

  return tempTagFile;
  }


TAGFILE *newTagFileFromFileState( FILESTATE *fileState )
  {
  TAGFILE *tempTagFile;

  tempTagFile = malloc( sizeof( TAGFILE ) );
  if( tempTagFile == NULL )
    {
    fprintf( stderr, "Memory allocation failure in "
      "newTagFileFromFileState()\n" );
    return NULL;
    }

  tempTagFile->currentTag = NULL;
  tempTagFile->file = fileState;
  tempTagFile->deallocateFileState = False;
  tempTagFile->readPending = True;

  return tempTagFile;
  }


void disposeTagFile( TAGFILE *tagFile )
  {
  if( tagFile->deallocateFileState == True )
    disposeFileState( tagFile->file );

  if( tagFile->currentTag != NULL )
    disposeTag( tagFile->currentTag );

  free( tagFile );
  }


void printTag( TAG *theTag )
  {
  if( theTag != NULL )
    printf( "%s: %s", theTag->tag, theTag->body );
  else 
    printf( "printTag(): NULL tag\n" );
  }


BOOLEAN tagCompare( const char *firstTag, const char *secondTag )
  {
  const char *s;
  const char *t;

  if( tagWhitespaceTableInitialized == 0 )
    {
    initWhitespaceTable( tagWhitespaceTable );
    tagWhitespaceTableInitialized = 1;
    }

  s = firstTag;
  t = secondTag;

  while( *s != '\0' && *t != '\0' )
    {
    while( tagWhitespaceTable[ (unsigned int) *s ] == TagWhitespace )
      s++;
    while( tagWhitespaceTable[ (unsigned int) *t ] == TagWhitespace )
      t++;

    if( tolower( *s++ ) != tolower ( *t++ ) )
      return False;
    }

  while( tagWhitespaceTable[ (unsigned int) *s ] == TagWhitespace )
      s++;
  while( tagWhitespaceTable[ (unsigned int) *t ] == TagWhitespace )
      t++;

  if( *s != '\0' || *t != '\0' )
    return False;
  else
    return True;
  }


FILEPOS getTagFilePosition( TAGFILE *theFile )
  {
  if( theFile->readPending == True )
    (void) getTag( theFile );

  return theFile->currentFilePos;
  }


void setTagFilePosition( TAGFILE *theFile, FILEPOS thePosition )
  { 
  setFilePosition( theFile->file, thePosition );
  theFile->currentFilePos = thePosition;
  theFile->readPending = True;
  } 
