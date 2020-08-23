/*
 *  Collection.h
 *
 *  The Collection module - Public Interface
 */

#ifndef COLLECTION_H
#define COLLECTION_H


/* -- Required Include Files --------------------------------------------------
*/

#include "groups.h"
    /* For: GROUPLIST below */

#include "reducers.h"
    /* For: reduce*() */



/* -- Public Module Datatypes -------------------------------------------------
*/

typedef struct COLLECTIONGROUPLIST
  {
  char                   *type;
  GROUPLIST              *groupList;

  struct COLLECTIONGROUPLIST *nextCollectionGroupList;
  }
COLLECTIONGROUPLIST;


typedef struct COLLECTION
  {
  COLLECTIONGROUPLIST *firstCollectionGroupList;
  }
COLLECTION;



/* -- Constructors and Destructors --------------------------------------------
*/ 

COLLECTIONGROUPLIST *newCollectionGroupList( char *type );
void disposeCollectionGroupList( COLLECTIONGROUPLIST *collectionGroupList );
void disposeCollectionGroupLists( COLLECTIONGROUPLIST *collectionGroupList );
COLLECTION *newCollection( void );
void disposeCollection( COLLECTION *collection );



/* -- Public Methods ----------------------------------------------------------
*/


GROUPLIST *lookupGroupList( COLLECTION *collection, char *type );
void addGroupToCollection( GROUP *group, COLLECTION *collection );
COLLECTION *newCollectionFromFile( char *baseName );
void writeCollection( COLLECTION *collection, char *baseName );
GROUP *lookupReducedGroup( GROUPLIST *groupList, char *name, 
  REDUCTION_FUNCTION *reduce );


#endif
