/*
 *  Collection.c
 *
 *  The Collection module - a container for lists of groups.
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stderr */

#include <errno.h>

#include <stdlib.h>
    /* For: malloc(), free(), exit(), NULL */

#include <string.h>
    /* For: strlen() */

#include "lineio.h"
    /* For: FILESTATE */

#include "tags.h"
    /* For: TAG, getTag() */

#include "groups.h"
    /* For: GROUP */

#include "rules.h"
    /* For: RULES, newRules(), etc */

#include "rule-parse.h"
    /* For: parseRuleFile() etc */

#include "utilities.h"
    /* For: outputString() */

#include "reducers.h"
    /* For: reduce*() */

#include "collection.h" 
    /* For: Our public interface */



/* -- Constructors and Destructors --------------------------------------------
*/ 

COLLECTIONGROUPLIST *newCollectionGroupList( char *type )
  {
  COLLECTIONGROUPLIST *tempCollectionGroupList;

  tempCollectionGroupList = malloc( sizeof( COLLECTIONGROUPLIST ) );
  if( tempCollectionGroupList == NULL )
    {
    fprintf( stderr, "Memory allocation failure in "
      "newCollectionGroupList() [1]\n" );
    return NULL;
    }

  tempCollectionGroupList->type = malloc( strlen( type ) + 1 );
  if( tempCollectionGroupList->type == NULL )
    {
    fprintf( stderr, "Memory allocation failure in "
      "newCollectionGroupList() [2]\n" );
    free( tempCollectionGroupList );
    return NULL;
    }

  strcpy( tempCollectionGroupList->type, type );

  tempCollectionGroupList->groupList = newGroupList(); 
  if( tempCollectionGroupList->groupList == NULL )
    {
    fprintf( stderr, "newGroupList() returned NULL in "
      "newCollectionGroupList()\n" );
    free( tempCollectionGroupList->type );
    free( tempCollectionGroupList );
    return NULL;
    }
 
  tempCollectionGroupList->nextCollectionGroupList = NULL;

  return tempCollectionGroupList;
  }


void disposeCollectionGroupList( COLLECTIONGROUPLIST *collectionGroupList )
  {
  free( collectionGroupList->type );
  disposeGroupList( collectionGroupList->groupList );
  free( collectionGroupList );
  }


void disposeCollectionGroupLists( COLLECTIONGROUPLIST *collectionGroupList )
  {
  if( collectionGroupList != NULL )
    {
    disposeCollectionGroupLists( collectionGroupList->nextCollectionGroupList );
    disposeCollectionGroupList( collectionGroupList );
    }
  }


COLLECTION *newCollection( void )
  {
  COLLECTION *tempCollection;

  tempCollection = malloc( sizeof( COLLECTION ) );
  if( tempCollection == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newCollection()\n" );
    return NULL;
    }

  tempCollection->firstCollectionGroupList = NULL;

  return tempCollection;
  }


void disposeCollection( COLLECTION *collection )
  {
  if( collection->firstCollectionGroupList != NULL )
    disposeCollectionGroupLists( collection->firstCollectionGroupList );

  free( collection );
  }








/* -- Public Methods ----------------------------------------------------------
*/


GROUPLIST *lookupGroupList( COLLECTION *collection, char *type )
  {
  COLLECTIONGROUPLIST *collectionGroupList;

  collectionGroupList = collection->firstCollectionGroupList;
  while( collectionGroupList != NULL )
    {
    if( tagCompare( type, collectionGroupList->type ) == True )
      return collectionGroupList->groupList;
    
    collectionGroupList = collectionGroupList->nextCollectionGroupList;
    }

  return NULL;
  }


void addGroupToCollection( GROUP *group, COLLECTION *collection )
  { 
  GROUPLIST           *groupList;
  GROUP               *groupCopy;
  COLLECTIONGROUPLIST *collectionGroupList;
    
  groupList = lookupGroupList( collection, group->type );
  if( groupList == NULL )
    {
    collectionGroupList = newCollectionGroupList( group->type );
    if( collectionGroupList == NULL )
      {
      fprintf( stderr, "newOutputGroupList() returned NULL in "
        "addGroupToOutput()\n" );
      return;
      }
 
    collectionGroupList->nextCollectionGroupList =
      collection->firstCollectionGroupList;
    collection->firstCollectionGroupList = collectionGroupList;
    }

  collectionGroupList = collection->firstCollectionGroupList;
  while( collectionGroupList != NULL )
    {
    if( tagCompare( group->type, collectionGroupList->type ) == True )
      {
      groupCopy = duplicateGroup( group );

      if( groupCopy == NULL )
        {
        fprintf( stderr, "duplicateGroup() returned NULL in "
          "addGroupToOutput()\n" );
        return;
        }

      addGroupToGroupList( groupCopy, collectionGroupList->groupList );
      return;
      }

    collectionGroupList = collectionGroupList->nextCollectionGroupList;
    }
  }


COLLECTION *newCollectionFromFile( char *baseName )
  {
  COLLECTION *tempCollection;
  GROUPFILE  *groupFile;
  RULES      *rules;
  GROUP      *group;
  char       *tempFileName;

  tempCollection = newCollection();
  if( tempCollection == NULL )
    {
    fprintf( stderr, "newCollection() returned NULL in "
      "newCollectionFromFile()\n" );
    return NULL;
    }

  tempFileName = malloc( strlen( baseName ) + 1 + 6 );
  if( tempFileName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in readCollection()\n" );
    return NULL;
    }

  sprintf( tempFileName, "%s.list", baseName );
  groupFile = newGroupFile( tempFileName );
  if( groupFile == NULL )
    {
    fprintf( stderr, "newGroupFile() returned NULL in "
      "newCollectionFromFile()\n" );
    free( tempFileName );
    return NULL;
    }

  rules = newRules();
  if( rules == NULL )
    {
    fprintf( stderr, "newRules() returned NULL in "
      "newCollectionFromFile()\n" );
    disposeGroupFile( groupFile );
    free( tempFileName );
    return NULL;
    }

  sprintf( tempFileName, "%s.rules", baseName );
  parseRuleFile( rules, tempFileName );
    /* -- No error checking is preformed for this */

  group = getGroup( groupFile, rules );
  while( group != NULL )
    {
    addGroupToCollection( group, tempCollection );

    discardGroup( groupFile );
    group = getGroup( groupFile, rules );
    }

  return tempCollection;
  }


void outputTag( FILE *file, TAG *tag )
  {
  int numLines;
  int pendingNewlines;

  char *p;  /* pointer */
  char *s;  /* start */
  char *l;  /* last */

  /* Step 1 - Find the start */
  p = tag->body;
  while( *p != '\0' )
    {
    if( *p != '\n' && *p != '\r' && *p != '\t' && *p != ' ' )
      break;
    p++;
    }
  s = p;  /* s either equals the first non-whitespace char or end: '\0' */

  /* Step 2 - Find the end */
  l = p;  /* Set initial value */
  pendingNewlines = 0;
  numLines = 1;
  while( *p != '\0' )
    {
    if( *p == '\n' )
      {
      pendingNewlines++;
      }
    else
      {
      l = p;
      if( pendingNewlines != 0 )
        {
        numLines += pendingNewlines;
        pendingNewlines = 0;
        }
      }
    p++;
    }

  if( *l != '\0' )
    l++;   /* If l points to the last char and not the terminator,
              point it one more along (to the "sentinal" position */

  /* Step 3a - Output Tag */
  fprintf( file, "%s:", tag->tag );

  /* Step 3b - Output the above bounded string */
  if( numLines > 1 )
    {
    p = s;
    fprintf( file, "\n  " );
    while( p != l )
      {
      fputc( *p, file );
      if( *p == '\n' )
        fprintf( file, "  " );
      
      p++;
      }
    fprintf( file, "\n" );
    }
  else
    {
    fprintf( file, " " );
    fwrite( s, 1, l-s , file );
    fprintf( file, "\n" );
    }
  }


void outputGroup( FILE *file, GROUP *group )
  {
  TAG     *tag;
  TAGITEM *tagItem;

  tag = group->headTag;
  if( tag != NULL )
    outputTag( file, tag );

  tagItem = group->tagList->firstTagItem;
  while( tagItem != NULL )
    {
    tag = tagItem->tag;

    outputTag( file, tag );
    
    tagItem = tagItem->nextTagItem;
    }

  fprintf( file, "\n\n" );  /* Group Seperator */
  }


void outputGroupList( FILE *file, GROUPLIST *groupList )
  {
  GROUPITEM *groupItem;
  GROUP     *group;

  groupItem = groupList->firstGroupItem;
  while( groupItem != NULL )
    {
    group = groupItem->group;
    outputGroup( file, group );

    groupItem = groupItem->nextGroupItem;
    }
  }


void writeCollection( COLLECTION *collection, char *baseName )
  {
  FILE      *file;
  char      *fileName;
  GROUPLIST *groupList;
  COLLECTIONGROUPLIST *collectionGroupList;

  fileName = malloc( strlen( baseName ) + 7 );
  if( fileName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in writeCollection()\n" );
    exit( 1 );
    }

  sprintf( fileName, "%s.list", baseName );
  file = fopen( fileName, "wt" );
  if( file == NULL )
    {
    fprintf( stderr, "Failed to open file %s in writeCollection(): %s\n",
	     fileName, strerror(errno) );
    free( fileName );
    exit( 1 );
    }

  collectionGroupList = collection->firstCollectionGroupList;
  while( collectionGroupList != NULL )
    {
    groupList = collectionGroupList->groupList;
    outputGroupList( file, groupList );

    collectionGroupList = collectionGroupList->nextCollectionGroupList;
    }
  }


GROUP *lookupReducedGroup( GROUPLIST *groupList, char *name, 
  REDUCTION_FUNCTION *reduce )
  {
  GROUP     *group;
  GROUPITEM *groupItem;
  char       key[256];       /* Whoop whoop! Overflow about to happen! */
  char       groupName[256]; /*      ''          ''          ''        */

  reduce( key, name );

  groupItem = groupList->firstGroupItem;
  while( groupItem != NULL )
    {
    group = groupItem->group;
    reduce( groupName, group->name );
    if( strcmp( key, groupName ) == 0 )
      return group;
    
    groupItem = groupItem->nextGroupItem;
    }

  return NULL;
  }


TAG *newLineTag( int start, int end )
  {
  TAG  *tempTag;
  char *tempString;
  
  tempTag = newTag();
  if( tempTag == NULL )
    {
    fprintf( stderr, "newTag() returned NULL in newLineTag()\n" );
    return NULL;
    }

  tempString = malloc( strlen( "Line" ) + 1 );
  if( tempString == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newLineTag() [1]\n" );
    disposeTag( tempTag );
    return NULL;
    }

  strcpy( tempString, "Line" );
  tempTag->tag = tempString;

  tempString = malloc( 18 );
  if( tempString == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newLineTag() [2]\n" );
    free( tempTag->tag );
    disposeTag( tempTag );
    return NULL;
    }

  sprintf( tempString, "%d-%d", start, end );
  tempTag->body = tempString;

  return tempTag;
  }



/* -- Module test code --------------------------------------------------------
*/

#ifdef TESTING

#define TestInput "/usr/httpd/www.esc.auckland.ac.nz/Groups/Bioengineering/CMISS/help/documentation/code/variables"

int main( void )
  {
  COLLECTION *collection;
  GROUPLIST  *groupList;
  GROUP      *group;
  TAG        *lineTag;

  collection = newCollectionFromFile( TestInput );
  if( collection == NULL )
    {
    fprintf( stderr, "newCollectionFromFile() returned NULL in main()\n" );
    exit( 1 );
    }

  groupList = lookupGroupList( collection, "VARIABLE" );
  if( groupList == NULL )
    {
    fprintf( stderr, "lookupGroupList() returned NULL in main()\n" );
    exit( 1 );
    }

  group = lookupReducedGroup( groupList, "YD", reduceVariable );
  if( group == NULL )
    {
    fprintf( stderr, "lookupGroup() returned NULL in main()\n" );
    exit( 1 );
    }

  lineTag = newLineTag( 1000, 2000 );
  if( lineTag == NULL )
    {
    fprintf( stderr, "newTag() returned NULL in main()\n" );
    exit( 1 );
    }

  addTagToGroup( lineTag, group );

  writeCollection( collection, "blort" );

  disposeCollection( collection );

  return 0;
  }
#endif
