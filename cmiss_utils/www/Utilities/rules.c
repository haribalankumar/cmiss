/*
 *  Rules.c
 *
 *    Rules module
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), etderr, etc */

#include <stdlib.h>
    /* For: malloc(), free(), NULL */

#include <string.h>
    /* For: strcmp() */

#include "rules.h"
    /* For: our public interface */



/* -- Constructors and Destructors --------------------------------------------
*/

TAGRULE *newTagRule( char *name, int minOccurances, int maxOccurances )
  {
  TAGRULE *tempTagRule;

  tempTagRule = malloc( sizeof( TAGRULE ) );
  if( tempTagRule == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newTagRule()\n" );
    return NULL;
    }
  
  tempTagRule->headRule = False;
  tempTagRule->name = name;
  tempTagRule->minOccurances = minOccurances;
  tempTagRule->maxOccurances = maxOccurances;
  tempTagRule->nextTagRule = NULL;

  return tempTagRule;
  }


TAGRULE *newHeadTagRule( char *name, int minOccurances, int maxOccurances )
  {
  TAGRULE *tempTagRule;

  tempTagRule = newTagRule( name, minOccurances, maxOccurances );
  if( tempTagRule == NULL )
    {
    fprintf( stderr, "newTagRule() returned NULL in newHeadTagRule()\n" );
    return NULL;
    }

  tempTagRule->headRule = True;

  return tempTagRule;
  }


void disposeTagRule( TAGRULE *theTagRule )
  {
  free( theTagRule->name );
  free( theTagRule );
  }


void disposeTagRuleList( TAGRULE *theTagRule )
  {
  if( theTagRule != NULL )
    {
    disposeTagRuleList( theTagRule->nextTagRule );
    disposeTagRule( theTagRule );
    }
  }


GROUPRULE *newGroupRule( char *groupName )
  {
  GROUPRULE *tempGroupRule;

  tempGroupRule = malloc( sizeof( GROUPRULE ) );
  if( tempGroupRule == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newGroupRule()\n" );
    return NULL;
    }
  
  tempGroupRule->name = groupName;
  tempGroupRule->firstTagRule = NULL;
  tempGroupRule->headTagRule = NULL;
  tempGroupRule->nextGroupRule = NULL;

  return tempGroupRule;
  }


void disposeGroupRule( GROUPRULE *theGroupRule )
  {
  free( theGroupRule->name );
  disposeTagRule( theGroupRule->headTagRule );
  disposeTagRuleList( theGroupRule->firstTagRule );
  free( theGroupRule );
  }


RULES *newRules( void )
  {
  RULES *tempRules;

  tempRules = malloc( sizeof ( RULES ) );
  if( tempRules == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newRules()\n" );
    return NULL;
    }

  tempRules->firstGroupRule = NULL;

  return tempRules;
  }


void disposeRules( RULES *rules )
  {
  if( rules != NULL )
    free( rules );
  }



/* -- Query Methods -----------------------------------------------------------
*/

BOOLEAN headTag( TAG *tag, RULES *rules )
  {
  GROUPRULE *group;

  group = rules->firstGroupRule;
  while( group != NULL )
    {
    if( tagCompare( group->name, tag->tag ) == True )
      return True;

    group = group->nextGroupRule;
    }

  /* Failed to find a match */
  return False;
  }



/* -- Public Methods ----------------------------------------------------------
*/

void addGroupRuleToRules( RULES *rules, char *groupName )
  {
  GROUPRULE *groupRule;

  /* Should check to make sure an item doesn't already exist */

  groupRule = newGroupRule( groupName );
  if( groupRule == NULL )
    {
    fprintf( stderr, "newGroupRule() return NULL in addGroupToRules()\n" );
    return;
    }  
  
  groupRule->nextGroupRule = rules->firstGroupRule;
  rules->firstGroupRule = groupRule;
  }


void addHeadToGroupRule( RULES *rules, char *groupName, char *headName, 
  int minOccurances, int maxOccurances )
  {
  GROUPRULE *groupRule;
  TAGRULE   *headRule;

  groupRule = lookupGroupRule( rules, groupName );
  if( groupRule == NULL )
    {
    fprintf( stderr, "Group Rule does not exist in addHeadToGroupRule()\n" );
    return;
    }
  
  headRule = newHeadTagRule( headName, minOccurances, maxOccurances );
  if( headRule == NULL )
    {
    fprintf( stderr, "newHeadTagRule() returned NULL in "
      "addHeadToGroupRule()\n" );
    return;
    }
    
  /* We should check there isn't one already */

  groupRule->headTagRule = headRule;
  }
  

void addTagToGroupRule( RULES *rules, char *groupName, char *tagName, 
  int minOccurances, int maxOccurances )
  {
  GROUPRULE *groupRule;
  TAGRULE   *tagRule;

  groupRule = lookupGroupRule( rules, groupName );
  if( groupRule == NULL )
    {
    fprintf( stderr, "Group Rule does not exist in addTagToGroupRule()\n" );
    return;
    }

  tagRule = newTagRule( tagName, minOccurances, maxOccurances );
  if( tagRule == NULL )
    {
    fprintf( stderr, "newTagRule() returned NULL in addHeadToGroup()\n" );
    return;
    }
    
  /* We should check there isn't one by the same name already */

  tagRule->nextTagRule = groupRule->firstTagRule;
  groupRule->firstTagRule = tagRule;
  }



/* -- Debugging Methods -------------------------------------------------------
*/

void printRules( RULES *rules )
  {
  GROUPRULE *group;
  TAGRULE   *tag;

  group = rules->firstGroupRule;
  while( group != NULL )
    {
    printf( "define %s as\n  <\n", group->name );
    printf( "  HEAD %s\n", group->headTagRule->name );
    tag = group->firstTagRule;
    while( tag != NULL )
      {
      printf( "  TAG  %s", tag->name );
      if( tag->minOccurances != 1 || tag->maxOccurances != 1 )
        {
        printf( "[" );
        if( tag->minOccurances == UnboundedOccurances )
          printf( "N" );
        else
          printf( "%d", tag->minOccurances );
        printf( ".." );
        if( tag->maxOccurances == UnboundedOccurances )
          printf( "N" );
        else
          printf( "%d", tag->maxOccurances );
        printf( "]\n" );
        }
      else
        {
        printf( "\n" );
        }
      tag = tag->nextTagRule;
      }
    printf( "  >\n" );
    printf( "\n" );
    group = group->nextGroupRule;
    }
  }



/* -- Database (lookup) Methods -----------------------------------------------
*/

GROUPRULE *findGroupRuleByHead( TAG *tag, RULES *rules )
  {
  GROUPRULE *group;

  group = rules->firstGroupRule;

  while( group != NULL )
    {
    if( tagCompare( group->headTagRule->name, tag->tag ) == True )
      return group;

    group = group->nextGroupRule;
    }
    
  /* Failed to find match */
    
  return NULL;
  }


TAGRULE *findTagRuleInGroupRule( TAG *tag, GROUPRULE *group )
  {
  TAGRULE *rule;

  rule = group->firstTagRule;
  while( rule != NULL )
    {
    if( tagCompare( tag->tag, rule->name ) == True )
      return rule;

    rule = rule->nextTagRule;
    }

  /* Failed to find a match */
  return NULL;
  }


GROUPRULE *lookupGroupRule( RULES *rules, char *groupName )
  {
  GROUPRULE *groupRule;

  groupRule = rules->firstGroupRule;
  while( groupRule != NULL )
    {
    if( strcmp( groupName, groupRule->name ) == 0 )
      return groupRule;

    groupRule = groupRule->nextGroupRule;
    }

  return NULL;  /* Failed to find a match */
  }


