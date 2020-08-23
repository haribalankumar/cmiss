/*
 *  Reducers.c
 *
 *    String reduction methods (one ofr variables, one for commands, etc)
 */


/* -- Include Directives ------------------------------------------------------
*/

#include <ctype.h>
    /* For: isprint(), tolower() */
#include <string.h>
    /* For: isprint(), tolower() */

#include "reducers.h"
    /* For: our public interface */



/* -- Module Public Functions -------------------------------------------------
*/


/* 
 *  reduceVariable:
 *
 *  +  Removes Whitespace
 *  +  Stops at first '('
 *  +  Lowercases output
 *
 */

void reduceVariable( char *destination, const char *source )
  {
  char       *d;
  const char *s;

  d = destination;
  s = source;

  while( *s != '\0' )
    {
    if( *s == '(' )
      break;

    if( *s == ' ' )
      {
      s++;
      continue;
      }

    if( isprint( *s ) )
      {
      *d = tolower( *s );
      d++;
      }
      
    s++;
    }

  *d = '\0';
  }



/*
 *  reduceModule:
 *
 *  +  Removes Whitespace
 *  +  Lowercases output
 *  +  Removes .f and .cmn extensions
 */

void reduceModule( char *destination, const char *source )
  {
  char       *d, *dot;
  const char *s;
  
  d = destination;
  s = source;
  
  while( *s != '\0' )
    {
    if( *s == ' ' )
      {
      s++;
      continue;
      }
    if( isprint( *s ) )
      {
      *d = tolower( *s );
      d++;
      }
      
    s++;
    }
  
  *d = '\0';

  dot = strrchr( destination, '.' );
  if( dot != NULL
      && ( 0 == strcmp( dot, ".f" ) ||
	   0 == strcmp( dot, ".cmn" ) ) )
    *dot = '\0';
  }


/* 
 *  reduceBlockData:
 *
 *  +  Removes Whitespace
 *  +  Lowercases output
 *
 */

void reduceBlockData( char *destination, const char *source )
  {
  char       *d;
  const char *s;

  d = destination;
  s = source;

  while( *s != '\0' )
    {
    if( *s == ' ' )
      {
      s++;
      continue;
      }

    if( isprint( *s ) )
      {
      *d = tolower( *s );
      d++;
      }

    s++;
    }

  *d = '\0';
  }


/* 
 *  reduceSubroutine:
 *
 *  +  Removes Whitespace
 *  +  Lowercases output
 *
 */

void reduceSubroutine( char *destination, const char *source )
  {
  char       *d;
  const char *s;

  d = destination;
  s = source;

  while( *s != '\0' )
    {
    if( *s == ' ' )
      {
      s++;
      continue;
      }

    if( isprint( *s ) )
      {
      *d = tolower( *s );
      d++;
      }

    s++;
    }

  *d = '\0';
  }


/*
 *  reduceFunction:
 *
 *  +  Removes Whitespace
 *  +  Lowercases output
 *
 */

void reduceFunction( char *destination, const char *source )
  {
  char       *d;
  const char *s;

  d = destination;
  s = source;

  while( *s != '\0' )
    {
    if( *s == ' ' )
      {
      s++;
      continue;
      }

    if( isprint( *s ) )
      {
      *d = tolower( *s );
      d++;
      }

    s++;
    }

  *d = '\0';
  }



/*
 *  reduceComment:
 *
 *  +  Removes Whitespace
 *  +  Lowercases output
 *
 */

void reduceComment( char *destination, const char *source )
  {
  char       *d;
  const char *s;

  d = destination;
  s = source;

  while( *s != '\0' )
    {
    if( *s == ' ' )
      {
      s++;
      continue;
      }

    if( isprint( *s ) )
      {
      *d = tolower( *s );
      d++;
      }

    s++;
    }

  *d = '\0';
  }




/*
 *  reduceCommand:
 *
 *  +  stop at '/'
 *  +  remove patterns matching "<.*>"
 *  +  remove patterns matching "[.*]" 
 *  +  remove patterns matching ";[a-z/]*
 *  +  remove space & '-'
 *  +  stop at any all uppercase word
 *  +  stop at any word it feels like
 *
 *  Basically a real mess - the command syntax is human readable and
 *  human parsable (i.e. each command has code that parses it, written by 
 *  hand) instead of machine readable. Much in need of rationalisation
 *  via a formal grammar. (yacc would be nice, but only because it's less
 *  work).
 *
 *  None of this is implemented at the moment of course - I'll just 
 *  remove the spaces and give up at the first [;<[/] for now
 */

void reduceCommand( char *destination, const char *source )
  {
  char       *d;
  const char *s;

  d = destination;
  s = source;

  while( *s != '\0' )
    {
    if( *s == ' ' || *s == '-' )
      {
      s++;
      continue;
      }

/*      if( *s == '[' || *s == ';' || *s == '<' || *s == '/' )  */
/*        break;  */
      if( *s == '[' )
        break;  

    if( isprint( *s ) )
      {
      *d = tolower( *s );
      d++;
      }

    s++;
    }

  *d = '\0';
  }



