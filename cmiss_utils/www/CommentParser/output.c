/*
 *  Output.c
 *
 *    Parseing output module
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stderr */

#include <stdlib.h>
    /* For: NULL, malloc(), free() */

#include <string.h>
    /* For: strlen() */

#include "output.h"
    /* For: our public interface */

#include "groups.h"
    /* For: GROUP, GROUPLIST, newGroupList(), etc */



/* -- Public Constructors and Destructors -------------------------------------
*/

OUTPUTGROUPLIST *newOutputGroupList( char *type )
  {
  OUTPUTGROUPLIST *tempOutputGroupList;

  tempOutputGroupList = malloc( sizeof( OUTPUTGROUPLIST ) );
  if( tempOutputGroupList == NULL )
    {
    fprintf( stderr, "Memory allocation failure in "
      "newOutputGroupList() [1]\n" );
    return NULL;
    }

  tempOutputGroupList->type = malloc( strlen( type ) + 1 );
  if( tempOutputGroupList->type == NULL )
    {
    fprintf( stderr, "Memory allocation failure in "
      "newOutputGroupList() [2]\n" );
    free( tempOutputGroupList );
    return NULL;
    }

  tempOutputGroupList->groupList = newGroupList();
  if( tempOutputGroupList->groupList == NULL )
    {
    fprintf( stderr, "newGroupList() returned NULL in newOutputGroupList()\n" );
    free( tempOutputGroupList->type );
    free( tempOutputGroupList );
    return NULL;
    }

  strcpy( tempOutputGroupList->type, type );

  tempOutputGroupList->next = NULL;

  return tempOutputGroupList;
  }


void disposeOutputGroupList( OUTPUTGROUPLIST *outputGroupList )
  {
  disposeGroupList( outputGroupList->groupList );
  free( outputGroupList->type );
  free( outputGroupList );
  }


OUTPUT *newOutput( void )
  {
  OUTPUT *tempOutput;

  tempOutput = malloc( sizeof( OUTPUT ) );
  if( tempOutput == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newOutput()\n" );
    return NULL;
    }

  tempOutput->firstOutputGroupList = NULL;

  return tempOutput;
  }


void disposeOutput( OUTPUT *output )
  {
  if( output->firstOutputGroupList != NULL )
    disposeOutputGroupList( output->firstOutputGroupList );
               /* -- This doesn't free all it should?? Check and fix. */

  free( output );
  }



/* -- Public Methods ----------------------------------------------------------
*/

void addGroupToOutput( GROUP *group, OUTPUT *output )
  {
  GROUPLIST       *groupList;
  GROUP           *groupCopy;
  OUTPUTGROUPLIST *outputGroupList;
    
  groupList = lookupGroupList( output, group->type );
  if( groupList == NULL )
    {
    outputGroupList = newOutputGroupList( group->type );
    if( outputGroupList == NULL )
      {
      fprintf( stderr, "newOutputGroupList() returned NULL in "
        "addGroupToOutput()\n" );
      return;
      }
 
    outputGroupList->next = output->firstOutputGroupList;
    output->firstOutputGroupList = outputGroupList;
    }

  outputGroupList = output->firstOutputGroupList;
  while( outputGroupList != NULL )
    {
    if( tagCompare( group->type, outputGroupList->type ) == True )
      {
      groupCopy = duplicateGroup( group );

      if( groupCopy == NULL )
        {
        fprintf( stderr, "duplicateGroup() returned NULL in "
          "addGroupToOutput()\n" );
        return;
        }

      addGroupToGroupList( groupCopy, outputGroupList->groupList );
      return;
      }

    outputGroupList = outputGroupList->next;
    }
  }


GROUPLIST *lookupGroupList( OUTPUT *output, char *type )
  {
  OUTPUTGROUPLIST *outputGroupList;

  outputGroupList = output->firstOutputGroupList;
  while( outputGroupList != NULL )
    {
    if( tagCompare( type, outputGroupList->type ) == True )
      return outputGroupList->groupList;
    
    outputGroupList = outputGroupList->next;
    }

  return NULL;
  }
