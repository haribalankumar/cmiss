/*
 *  Outputcommands.h
 *
 *    Public interface to the Outputcommands module
 */

#ifndef OUTPUTCOMMANDS_H
#define OUTPUTCOMMANDS_H


/* -- Required Include Files --------------------------------------------------
*/

#include "groups.h" 
    /* For: GROUPS below */


typedef enum
  {
  CommaSeperatedList,
  UnorderedBulletList
  }
LISTFORMAT;


/* -- Public Datatypes --------------------------------------------------------
*/


/* -- Public Methods ----------------------------------------------------------
*/

/* -- Module Public Functions -------------------------------------------------
*/

void outputCommand( FILE *output, GROUP *group );
void outputCommand2( FILE *output, GROUP *group );
void outputCommand3( FILE *output, GROUP *group );
/*void outputNameList( FILE *output, char *string, LISTFORMAT format )  */

#endif
 
