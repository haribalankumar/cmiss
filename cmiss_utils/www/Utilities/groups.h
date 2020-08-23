/*
 *  Groups.h
 * 
 *    Public interface to the Groups module.
 */
#ifndef GROUPS_H
#define GROUPS_H


/* -- Required Includes -------------------------------------------------------
*/

#include "lineio.h"
   /* For: FILESTATE below */

#include "rules.h"
    /* For: RULES below */

#include "tags.h"
    /* For: TAG, TAGFILE, TAGLIST below */



/* -- Public Module Datatypes -------------------------------------------------
*/

typedef struct
  {
  char    *type;    /* The text name for the group from the rules */
  char    *name;    /* The body of the HEAD tag */
  TAG     *headTag; /* Head tag */
  TAGLIST *tagList;
  int      numTags;
  }
GROUP;


typedef struct GROUPITEM
  {
  GROUP *group;

  struct GROUPITEM *nextGroupItem;
  }
GROUPITEM;


typedef struct
  {
  int        numGroups;
  GROUPITEM *firstGroupItem;
  }
GROUPLIST;


typedef struct
  {
  GROUP   *currentGroup;
  FILEPOS  currentFilePos;     /* File position for get*FilePosition() */
  TAGFILE *tagFile;
  BOOLEAN  deallocateTagFile;  /* Do we manage the TAGFILE, or someone else? */
  BOOLEAN  deallocateFileState; /* Ditto for FILESTATE */
  BOOLEAN  readPending;
  }
GROUPFILE;



/* -- Public Constructors & Destructors ---------------------------------------
*/

GROUP *newGroup( char *groupType );
GROUP *newGroupFromHeadTag( TAG *headTag, char *groupType );
GROUP *duplicateGroup( GROUP *original );
void disposeGroup( GROUP *group );

GROUPITEM *newGroupItem( GROUP *group );
void disposeGroupItem( GROUPITEM *groupItem );

GROUPLIST *newGroupList( void );
void disposeGroupList( GROUPLIST *groupList );

GROUPFILE *newGroupFile( char *fileName );
GROUPFILE *newGroupFileFromFileState( FILESTATE *fileState );
GROUPFILE *newGroupFileFromTagFile( TAGFILE *tagFile );
void disposeGroupFile( GROUPFILE *groupFile );



/* -- Public Methods ----------------------------------------------------------
*/

void addHeadTagToGroup( TAG *headTag, GROUP *group );
void addTagToGroup( TAG *tag, GROUP *group );

GROUP *getGroup( GROUPFILE *groupFile, RULES *rules );
void discardGroup( GROUPFILE *groupFile );

FILEPOS getGroupFilePosition( GROUPFILE *groupFile );
void setGroupFilePosition( GROUPFILE *groupFile, FILEPOS thePosition );

void addGroupToGroupList( GROUP *group, GROUPLIST *groupList );

void printGroup( GROUP *group );

TAG *lookupTag( GROUP *group, char *type );
TAG *lookupNthTag( GROUP *group, char *type, int N );
int getNumberOfTags( GROUP *group, char *type );

GROUP *lookupGroup( GROUPLIST *groupList, char *name );

#endif
