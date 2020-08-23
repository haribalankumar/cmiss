/*
 *  Output.h
 *
 *    Public interface to the Parseing Output module
 */

#ifndef OUTPUT_H
#define OUTPUT_H


/* -- Required Include Files --------------------------------------------------
*/

#include "groups.h"
    /* For: GROUPS below */



/* -- Public Datatypes --------------------------------------------------------
*/

typedef struct OUTPUTGROUPLIST
  {
  char                   *type;
  GROUPLIST              *groupList;

  struct OUTPUTGROUPLIST *next;
  }
OUTPUTGROUPLIST;


typedef struct
  {
  OUTPUTGROUPLIST *firstOutputGroupList;
  }
OUTPUT;



/* -- Public Constructors and Destructors -------------------------------------
*/

OUTPUT *newOutput( void );
void disposeOutput( OUTPUT *output );



/* -- Public Methods ----------------------------------------------------------
*/

/* Group/Output manipulation */

void addGroupToOutput( GROUP *group, OUTPUT *output );
GROUP *lookupOutputGroup( char *type, char *key );
void sortOutputGroup( char *type );

/* Database type methods */

GROUPLIST *lookupGroupList( OUTPUT *output, char *type );


#endif
