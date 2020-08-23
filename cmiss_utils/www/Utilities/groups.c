/*
 *  Groups.c
 * 
 *    A Module for handling groups of tags.
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: printf(), etc */

#include <stdlib.h>
    /* For: malloc(), free(), NULL, etc */

#include <string.h>
    /* For: strcmp() */

#include "rules.h"
    /* For: RULES, newRules(), UnboundedOccurances etc */

#include "lineio.h"
    /* For: newFileState(), FILESTATE */

#include "tags.h"
    /* For: getTag() */

#include "groups.h"
    /* For: Our public interface */

#include "rule-parse.h"
    /* For: findGroupRuleByHead() */



/* -- Constructors and Destructors --------------------------------------------
*/

GROUP *newGroup( char *groupType )
  {
  GROUP *tempGroup;

  tempGroup = malloc( sizeof( GROUP ) );
  if( tempGroup == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newGroup() [1]\n" );
    return NULL;
    }

  tempGroup->type = malloc( strlen( groupType ) + 1 );
  if( tempGroup->type == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newGroup() [2]\n" );
    free( tempGroup );
    return NULL;
    }

  strcpy( tempGroup->type, groupType );
  tempGroup->name = NULL;
  tempGroup->headTag = NULL;
  tempGroup->tagList = newTagList();
  tempGroup->numTags = 0;

  return tempGroup;
  }


GROUP *newGroupFromHeadTag( TAG *head, char *groupType )
  {
  GROUP *tempGroup;

  tempGroup = newGroup( groupType );
  if( tempGroup == NULL )
    {
    fprintf( stderr, "newGroup() returned NULL in newGroupFromHeadTag()" );
    return NULL;
    }

  addHeadTagToGroup( head, tempGroup );

  return tempGroup;
  }


GROUP *duplicateGroup( GROUP *original )
  {
  GROUP *copy;
  TAG   *headCopy;

  headCopy = duplicateTag( original->headTag );
  if( headCopy == NULL )
    {
    fprintf( stderr, "duplicateTag() returned NULL in duplicateGroup()\n" );
    return NULL;
    }

  copy = newGroupFromHeadTag( headCopy, original->type );
  if( copy == NULL )
    {
    fprintf( stderr, "newGroup() returned NULL in duplicateGroup()\n" );
    return NULL;
    }

  copy->tagList = duplicateTagList( original->tagList );
  if( copy->tagList == NULL )
    {
    fprintf( stderr, "duplicateTagList() returned NULL in "
      "duplicateGroup()\n" );
    return NULL;
    }

  copy->numTags = original->numTags;

  return copy;
  }


void disposeGroup( GROUP *group )
  {
  if( group->name != NULL )
    free( group->name );

  if( group->headTag != NULL )
    disposeTag( group->headTag );
  
  disposeTagList( group->tagList );

  free( group );
  }


GROUPITEM *newGroupItem( GROUP *group )
  {
  GROUPITEM *tempGroupItem;

  tempGroupItem = malloc( sizeof( GROUPITEM ) );
  if( tempGroupItem == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newGroupItem()\n" );
    return NULL;
    }

  tempGroupItem->group = group;
  tempGroupItem->nextGroupItem = NULL;

  return tempGroupItem;
  }


void disposeGroupItem( GROUPITEM *groupItem )
  {
  disposeGroup( groupItem->group );
  free( groupItem ); 
  }


void disposeGroupItemList( GROUPITEM *groupItem )
  {
  if( groupItem != NULL )
    {
    disposeGroupItemList( groupItem->nextGroupItem );
    disposeGroupItem( groupItem );
    }
  }


GROUPLIST *newGroupList( void )
  {
  GROUPLIST *tempGroupList;

  tempGroupList = malloc( sizeof( GROUPLIST ) );
  if( tempGroupList == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newGroupList()\n" );
    return NULL;
    }

  tempGroupList->numGroups = 0;
  tempGroupList->firstGroupItem = NULL;

  return tempGroupList;  
  }


void disposeGroupList( GROUPLIST *groupList )
  {
  disposeGroupItemList( groupList->firstGroupItem );
  free( groupList );
  }


GROUPFILE *newGroupFile( char *fileName )
  {
  GROUPFILE *tempGroupFile;
  FILESTATE *fileState;

  fileState = newFileState( fileName );
  if( fileState == NULL )
    {
    fprintf( stderr, "newFileState() returned NULL in newGroupFile()\n" );
    return NULL;
    }

  tempGroupFile = newGroupFileFromFileState( fileState );
  if( tempGroupFile == NULL )
    {
    fprintf( stderr, "newGroupFileFromFileState() returned NULL in "
      "newGroupFile()\n" );
    disposeFileState( fileState );
    return NULL;
    }

  tempGroupFile->deallocateFileState = True;

  return tempGroupFile;
  }


GROUPFILE *newGroupFileFromFileState( FILESTATE *fileState )
  {
  GROUPFILE *tempGroupFile;
  TAGFILE   *tagFile;

  tagFile = newTagFileFromFileState( fileState );
  if( tagFile == NULL )
    {
    fprintf( stderr, "newTagFileFromFileState() returned NULL in "
      "newGroupFileFromFileState()\n" );
    return NULL;
    }

  tempGroupFile = newGroupFileFromTagFile( tagFile );
  if( tempGroupFile == NULL )
    {
    fprintf( stderr, "newGroupFileFromTagFile() returned NULL in "
      "newGroupFileFromFileState()\n ");
    disposeTagFile( tagFile );
    return NULL;
    }

  tempGroupFile->deallocateTagFile = True;  /* Override default value */

  return tempGroupFile;
  }


GROUPFILE *newGroupFileFromTagFile( TAGFILE *tagFile )
  {
  GROUPFILE *tempGroupFile;

  tempGroupFile = malloc( sizeof( GROUPFILE ) );
  if( tempGroupFile == NULL )
    {
    fprintf( stderr, "Memory allocation failure in "
      "newGroupFileFromTagFile()\n" );
    return NULL;
    }

  tempGroupFile->tagFile = tagFile;
  tempGroupFile->currentGroup = NULL;
  tempGroupFile->currentFilePos = getTagFilePosition( tempGroupFile->tagFile );
      /* -- Unnessesery - should be overwritten later */
  tempGroupFile->deallocateTagFile = False;
  tempGroupFile->deallocateFileState = False;
  tempGroupFile->readPending = True;

  return tempGroupFile;
  }


void disposeGroupFile( GROUPFILE *groupFile )
  {
  if( groupFile->deallocateFileState == True )
    if( groupFile->tagFile != NULL )
      if( groupFile->tagFile->file != NULL )
        disposeFileState( groupFile->tagFile->file );

  if( groupFile->deallocateTagFile == True )
    if( groupFile->tagFile != NULL )
      disposeTagFile( groupFile->tagFile );

  if( groupFile->currentGroup != NULL )
    disposeGroup( groupFile->currentGroup );
 
  free( groupFile );
  }



/* -- Methods -----------------------------------------------------------------
*/

void addHeadTagToGroup( TAG *headTag, GROUP *group )
  {
  group->name = malloc( strlen( headTag->body ) + 1 );
  strcpy( group->name, headTag->body );

  group->headTag = headTag;
  }


void addTagToGroup( TAG *tag, GROUP *currentGroup )
  {
  TAGITEM *tagItem;
  TAGLIST *tagList;
  /* We should be checking for the validity of this tag somewhere - although
     this would require passing in a RULES structure */
  
  tagItem = newTagItem( tag );
  if( tagItem == NULL )
    {
    fprintf( stderr, "newTagItem() returned NULL in addTagToGroup()\n" );
    return;
    }

  tagList = currentGroup->tagList;

  if( tagList->lastTagItem != NULL )
    {
    tagList->lastTagItem->nextTagItem = tagItem;
    tagList->lastTagItem = tagItem;
    }
  else
    {
    tagList->firstTagItem = tagList->lastTagItem = tagItem;
    }

  tagList->numEntries++;
  }



/* -- Module Public Methods ---------------------------------------------------
*/ 

GROUP *getGroup( GROUPFILE *groupFile, RULES *rules )
  {
  GROUPRULE *groupRule;
  TAGRULE   *tagRule;
  TAG       *tag;
  TAG       *tagCopy;


  if( groupFile->readPending == True )
    {
    if( groupFile->currentGroup != NULL )
      {
      disposeGroup( groupFile->currentGroup );
      groupFile->currentGroup = NULL;   /* So people know whether to free */
      }

    /* 
     *  This while() loop should be shot. Rewrite without the multiple
     *  discardTag()'s and esp. without any break's or continue's.
     */
    tag = getTag( groupFile->tagFile );

    while( tag != NULL )
      {
      tag = getTag( groupFile->tagFile );

      if( tag == NULL ) /* end of input */
        {
        /* Tidy this mess up later */
        break;
        }
    
      groupRule = findGroupRuleByHead( tag, rules );
      if( groupRule == NULL )
        {
        fprintf( stderr, "Tag (%s) encounted by getGroup() does not match "
          "any group rule (skipping)\n", tag->tag ); 
        fflush( stderr );
        discardTag( groupFile->tagFile );

        continue;
        }
      else
        {
        groupFile->currentFilePos = getTagFilePosition( groupFile->tagFile );

        tagCopy = duplicateTag( tag );
        if( tagCopy == NULL )
          {
          fprintf( stderr, "duplicateTag() returned NULL in getGroup()\n" );

          break;  /* Exit the while loop */
          }

        groupFile->currentGroup =
          newGroupFromHeadTag( tagCopy, groupRule->name );

        if( groupFile->currentGroup == NULL )
          {   
          fprintf( stderr, "newGroupFromHeadTag() returned NULL in "
            "getGroup()\n" );
          return NULL;
          }

        discardTag( groupFile->tagFile );
        break;
        }

      /* Never Reached */
      /* discardTag( groupFile->tagFile ); */
      }

    if( tag == NULL )
      {
      /* We have encountered some error, or have run out of input */
      return NULL;
      }

    tag = getTag( groupFile->tagFile );
    while( tag != NULL )
      {
      tagRule = findTagRuleInGroupRule( tag, groupRule );
      
      if( tagRule != NULL )
        {
        tagCopy = duplicateTag( tag );
        if( tagCopy == NULL )
          {
          fprintf( stderr, "duplicateTag() returned NULL in getGroup()\n" );
          break;  /* Exit the while loop */
          }
        
        addTagToGroup( tagCopy, groupFile->currentGroup );
        
        discardTag( groupFile->tagFile );
        }
      else
        {
        break;
        }

      tag = getTag( groupFile->tagFile );
      }

    groupFile->readPending = False;
    }

  return groupFile->currentGroup;
  }


void discardGroup( GROUPFILE *groupFile )
  {
  groupFile->readPending = True;
  }


void printGroup( GROUP *group )
  {
  TAGITEM *tagItem;

  printf( "Group type is %s\n", group->type );
  printf( "Group name is %s\n", group->name );
  printf( "  Head = " ); 
  printTag( group->headTag );
  
  tagItem = group->tagList->firstTagItem;
  while( tagItem != NULL )
    {
    printf( "  Tag = " );
    printTag( tagItem->tag );
    tagItem = tagItem->nextTagItem;
    }
  }


TAG *lookupTag( GROUP *group, char *type )
  {
  TAGITEM *tagItem;

  /* first check the head tag */
  if( tagCompare( type, group->headTag->tag ) == True )
    {
    return group->headTag;
    }

  /* then check the body tags */
  tagItem = group->tagList->firstTagItem;
  while( tagItem != NULL )
    {
    if( tagCompare( type, tagItem->tag->tag ) == True )
      {
      return tagItem->tag;
      }
    
    tagItem = tagItem->nextTagItem;
    }

  return NULL;
  }


TAG *lookupNthTag( GROUP *group, char *type, int N )
  { 
  TAGITEM *tagItem;
  int      number = 0;

  /* first check the head tag */
  if( tagCompare( type, group->headTag->tag ) == True )
    number++;

  if( number == N )
    return group->headTag;

  /* then check the body tags */
  tagItem = group->tagList->firstTagItem;
  while( tagItem != NULL )
    {
    if( tagCompare( type, tagItem->tag->tag ) == True )
      number++;

    if( number == N )
      return tagItem->tag;

    tagItem = tagItem->nextTagItem;
    }

  /* failed to find item */
  return NULL;
  }


int getNumberOfTags( GROUP *group, char *type )
  {
  TAGITEM *tagItem;
  int      number = 0;

  /* first check the head tag */
  if( tagCompare( type, group->headTag->tag ) == True )
    number++;

  /* then check the body tags */
  tagItem = group->tagList->firstTagItem;
  while( tagItem != NULL )
    {
    if( tagCompare( type, tagItem->tag->tag ) == True )
      number++;

    tagItem = tagItem->nextTagItem;
    }

  return number;
  }


GROUP *lookupGroup( GROUPLIST *groupList, char *name )
  {
  GROUP     *group;
  GROUPITEM *groupItem;

  groupItem = groupList->firstGroupItem;
  while( groupItem != NULL )
    {
    group = groupItem->group;

    if( strcmp( name, group->name ) == 0 )
      return group;
    
    groupItem = groupItem->nextGroupItem;
    }

  return NULL;
  }


void addGroupToGroupList( GROUP *group, GROUPLIST *groupList )
  {
  GROUPITEM *groupItem;

  groupItem = newGroupItem( group );
  if( groupItem == NULL )
    {
    fprintf( stderr, "newGroupItem() returned NULL in "
      "addGroupToGroupList()\n" );
    return;
    }

  groupItem->nextGroupItem = groupList->firstGroupItem;
  groupList->firstGroupItem = groupItem;
  }


FILEPOS getGroupFilePosition( GROUPFILE *groupFile )
  {
  if( groupFile->readPending == True )
    {
    printf( "Error: cannot obtain Group File Position after a discardGroup() "
      "call.\n  Am returning wrong value instead :-)\n" ); 
    }

  return groupFile->currentFilePos;
  }


void setGroupFilePosition( GROUPFILE *groupFile, FILEPOS thePosition )
  {
  setTagFilePosition( groupFile->tagFile, thePosition );
  groupFile->currentFilePos = thePosition;
  groupFile->readPending = True;
  } 
