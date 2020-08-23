/*
 *  Simple Parser.h
 *
 *    Public interface to the simple parsing module.
 */
#ifndef SIMPLE_PARSER_H
#define SIMPLE_PARSER_H


/* -- Public Datatypes --------------------------------------------------------
*/

#if !defined( BOOLEAN )
#define BOOLEAN int
#define True 1
#define False 0
#endif

typedef struct
  {
  char *stream;
  char *currentPtr;
  char *end;
  }
PARSESTATE;



/* -- Constructors and Destructors --------------------------------------------
*/

PARSESTATE *newParseState( char *theStream );
void disposeParseState( PARSESTATE *theParseState );



/* -- Module Method Prototypes ------------------------------------------------
*/

BOOLEAN prefix( PARSESTATE *theState, char *thePrefix );
BOOLEAN empty( PARSESTATE *theState );
void shiftChars( PARSESTATE *theParseState, int shiftAmount );
void skip( PARSESTATE *theState, char *skipList );
void skipWord( PARSESTATE *theState, char *theWord );
void skipUntil( PARSESTATE *theState, char *terminationList );
char *extract( PARSESTATE *theState, char *terminationList );


#endif
