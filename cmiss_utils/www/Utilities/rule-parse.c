/*
 *  Rule Parse.c
 * 
 *    A module for parsing rules files
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: printf(), etc */

#include <errno.h>

#include <stdlib.h>
    /* For: malloc(), free(), NULL, etc */

#include <string.h>
    /* For: strcmp() */

#include "simple-parser.h"
    /* For: PARSESTATE, newParseState(), prefix(), skip(), extract(), etc */

#include "rules.h"
    /* For: RULES, newRules(), UnboundedOccurances etc */

#include "rule-parse.h"
    /* For: Our public interface */



/* -- Private Module Constants ------------------------------------------------
*/ 

#define WhiteSpace          " \t\n\r"
#define EntryTerminator     "[ \t\n\r"
#define OccuranceSeperator  ". \t\n\r"
#define OccuranceTerminator "] \t\n\r"



/* -- Private Method Prototypes -----------------------------------------------
*/

static void skipWhitespace( PARSESTATE *state );
static char *loadFileAsString( const char *fileName );
static void parseHead( RULES *rules, char *groupName, PARSESTATE *state );
static void parseTag( RULES *rules, char *groupName, PARSESTATE *state );
static void parseNumber( RULES *rules, char *groupName, PARSESTATE *state );
static void parseGroupList( RULES *rules, char *groupName, PARSESTATE *state );
static void parseDefinition( RULES *rules, PARSESTATE *state );
static void parse( PARSESTATE *state, RULES *rules );



/* -- Private Utility Methods -------------------------------------------------
*/

/*
 *  Whitespace is defined as spaces, tabs, newlines, or comments.
 */
static void skipWhitespace( PARSESTATE *state )
  {
  int spaceCount = 0;  /* Termination condition */

  while( spaceCount != 2 )
    {
    if( prefix( state, "--" ) == True )
      {
      skipUntil( state, "\n" );
      spaceCount = 0;
      }
    else
      {
      skip( state, WhiteSpace );
      spaceCount++;
      }
    }
  }


static char *loadFileAsString( const char *fileName )
                                  /* -- do this with a mmap() in future? */
  {
  FILE *file;
  long  length;
  char *string;

  file = fopen( fileName, "rb" );
  if( file == NULL )
    {
    fprintf( stderr, "Failed to open file %s in loadFileAsString(): %s\n",
	     fileName, strerror(errno) );
    return NULL;
    }
 
  /* There must be a better way of doing this */
  fseek( file, 0, SEEK_END );
  length = ftell( file );

  string = malloc( length + 1 );
  if( string == NULL )
    {
    fprintf( stderr, "Memory allocation failure in loadFileAsString()\n" );
    return NULL;
    }
 
  fseek( file, 0, SEEK_SET );
  fread( string, 1, length, file );
  string[length] = '\0';

  fclose( file );

  return string;
  }



/* -- Private Module Methods --------------------------------------------------
*/

static void parseHead( RULES *rules, char *groupName, PARSESTATE *state )
  {
  char *headName  = NULL;
  char *minString = NULL;
  char *maxString = NULL;
  int   minOccurances;
  int   maxOccurances;

  skipWhitespace( state );
  headName = extract( state, EntryTerminator );
  skipUntil( state, EntryTerminator );
  skipWhitespace( state );
  
  if( prefix( state, "[" ) == True )
    {
    skipWord( state, "[" );
    skipWhitespace( state );
    minString = extract( state, OccuranceSeperator );
    skipWhitespace( state );
    if( prefix( state, ".." ) != True )
      {
      fprintf( stderr, "'..' expected but not found in parseHead() "
        "(bailing out)\n" );
      free( headName );
      free( minString );
      return;
      }
    skipWord( state, ".." );
    skipWhitespace( state );
    maxString = extract( state, OccuranceTerminator );
    skipWhitespace( state );
    if( prefix( state, "]" ) != True )
      {
      fprintf( stderr, "']' expected bvut not found in parseHead() "
        "(ignoring)\n" );
      }
    skipWord( state, "]" );
    skipWhitespace( state );
     
    minOccurances = atoi( minString );
    
    if( strcmp( maxString, "N" ) == 0 )
      {
      maxOccurances = UnboundedOccurances;
      }
    else
      {
      maxOccurances = atoi( maxString );
      }

    free( minString );
    free( maxString );
    }
  else
    {    
    minOccurances = 1;
    maxOccurances = 1;
    }

  addHeadToGroupRule( rules, groupName, headName, minOccurances,
    maxOccurances );
  }


static void parseTag( RULES *rules, char *groupName, PARSESTATE *state )
  {
  char *tagName;
  char *minString;
  char *maxString;
  int   minOccurances;
  int   maxOccurances;

  skipWhitespace( state );
  tagName = extract( state, EntryTerminator );
  skipUntil( state, EntryTerminator );
  skipWhitespace( state );
  
  if( prefix( state, "[" ) == True )
    {
    skipWord( state, "[" );
    skipWhitespace( state );
    minString = extract( state, OccuranceSeperator );
    skipWhitespace( state );
    if( prefix( state, ".." ) != True )
      {
      fprintf( stderr, "'..' expected bvut not found in parseTag() "
        "(bailing out)\n" );
      free( tagName );
      free( minString );
      return;
      }
    skipWord( state, ".." );
    skipWhitespace( state );
    maxString = extract( state, OccuranceTerminator );
    skipWhitespace( state );
    if( prefix( state, "]" ) != True )
      {
      fprintf( stderr, "']' expected bvut not found in parseTag() "
        "(ignoring)\n" );
      }
    skipWord( state, "]" );
    skipWhitespace( state );
     
    minOccurances = atoi( minString );
    
    if( strcmp( maxString, "N" ) == 0 )
      {
      maxOccurances = UnboundedOccurances;
      }
    else
      {
      maxOccurances = atoi( maxString );
      }

    free( minString );
    free( maxString );
    }
  else
    {    
    minOccurances = 1;
    maxOccurances = 1;
    }

  addTagToGroupRule( rules, groupName, tagName, minOccurances, maxOccurances );
  }


static void parseNumber( RULES *rules, char *groupName, PARSESTATE *state )
  {
  printf( "please write parseNumber()\n" );
  }


static void parseGroupList( RULES *rules, char *groupName, PARSESTATE *state )
  {
  BOOLEAN endOfList;

  endOfList = False;

  while( endOfList != True )
    {
    /* HEAD */
    if( prefix( state, "HEAD" ) == True )
      {
      skipWord( state, "HEAD" );
      skipWhitespace( state );
      parseHead( rules, groupName, state );
      }

    /* TAG */
    if( prefix( state, "TAG" ) == True )
      {
      skipWord( state, "TAG" );
      skipWhitespace( state );
      parseTag( rules, groupName, state );
      }

    /* NUMBER */
    if( prefix( state, "NUMBER" ) == True )
      {
      skipWord( state, "NUMBER" );
      skipWhitespace( state );
      parseNumber( rules, groupName, state );
      }

    /* End of list */
    if( prefix( state, ">" ) == True )
      {
      /* Leave it for our caller to remove the '>' from the stream */
      endOfList = True;
      }

    /* End of input */
    if( empty( state ) == True )
      {
      endOfList = True;
      }
    }
  }


static void parseDefinition( RULES *rules, PARSESTATE *state )
  {
  char *groupName;

  /* define */
  skipWhitespace( state );
  if( prefix( state, "define" ) != True )
    {
    fprintf( stderr, "'define' expected (but not found) in "
      "parseDefinition() -- bailing out\n" );
    return;
    }
  skipWord( state, "define" );

  /* define SOMETHING */
  skipWhitespace( state );
  groupName = extract( state, WhiteSpace );
  skipUntil( state, WhiteSpace );

  addGroupRuleToRules( rules, groupName );

  /* define SOMETHING as */
  skipWhitespace( state );
  if( prefix( state, "as" ) != True )
    {
    fprintf( stderr, "'as' expected (but not found) in "
      "parseDefinition() -- bailing out\n" );
    return;
    }
  skipWord( state, "as" );

  /* define SOMETHING as < */
  skipWhitespace( state );
  if( prefix( state, "<" )  != True )
    {
    fprintf( stderr, "'<' expected (but not found) in "
      "parseDefinition() -- bailing out\n" );
    return;
    }
  skipWord( state, "<" );

  /* < ... > */
  skipWhitespace( state );
  parseGroupList( rules, groupName, state ); /* The actual work */
  skipWhitespace( state );
  if( prefix( state, ">" )  != True )
    {
    fprintf( stderr, "'>' expected (but not found) in "
      "parseDefinition() -- bailing out\n" );
    return;
    }
  skipWord( state, ">" );
  }

 
static void parse( PARSESTATE *state, RULES *rules )
  {
  while( empty( state ) != True )
    {
    parseDefinition( rules, state );
    skipWhitespace( state );
    }
  }



/* -- Public Module Methods ---------------------------------------------------
*/

int parseRuleFile( RULES *rules, const char *ruleFileName )
  {
  char       *file;  /* The file as a string in memory */
  PARSESTATE *state;

  file = loadFileAsString( ruleFileName );
  if( file == NULL )
    {
    fprintf( stderr, "loadFileAsString() returned NULL in parseRuleFile()\n" );
    return 0;
    }

  state = newParseState( file );
  if( state == NULL )
    {
    fprintf( stderr, "newParseState() return NULL in parseRuleFile()\n" );
    return 0;
    }

  parse( state, rules );

  return 1;
  }
