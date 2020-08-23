/*
 *  Rules.h
 *
 *    Rules modules public interface
 */

#ifndef RULES_H
#define RULES_H


/* -- Required Include Files --------------------------------------------------
*/

#include <limits.h>
    /* For: INT_MAX below */

#include "tags.h" 
    /* For: TAG (although I'm not sure I really want it here) */



/* -- Type definitions --------------------------------------------------------
*/

#if !defined( BOOLEAN )
#define BOOLEAN int
#define False 0
#define True 1
#endif


typedef struct TAGRULE
  {
  char    *name;
  int      minOccurances;
  int      maxOccurances;
  BOOLEAN  headRule;

  struct TAGRULE *nextTagRule;
  }
TAGRULE;


typedef struct GROUPRULE
  {
  char    *name;
  TAGRULE *headTagRule;
  TAGRULE *firstTagRule;

  struct GROUPRULE *nextGroupRule;
  }
GROUPRULE;


typedef struct
  {
  GROUPRULE *firstGroupRule;
  }
RULES;



/* -- Module Constants --------------------------------------------------------
*/ 

#define UnboundedOccurances INT_MAX



/* -- Module Ctors/Dtors ------------------------------------------------------
*/

RULES *newRules( void );
void disposeRules( RULES *rules );



/* -- Module Public Methods ---------------------------------------------------
*/

/* Adding rules */
void addGroupRuleToRules( RULES *rules, char *groupName );
void addHeadToGroupRule( RULES *rules, char *groupName, char* headName,
  int minOccurances, int maxOccurances );
void addTagToGroupRule( RULES *rules, char *groupName, char* tagName,
  int minOccurances, int maxOccurances );

/* Querying data */
BOOLEAN headTag( TAG *tag, RULES *rules );

/* Debugging */
void printRules( RULES *rules );

/* Database Lookup */
GROUPRULE *findGroupRuleByHead( TAG *tag, RULES *rules );
TAGRULE *findTagRuleInGroupRule( TAG *tag, GROUPRULE *group );
GROUPRULE *lookupGroupRule( RULES *rules, char *groupName );

#endif
