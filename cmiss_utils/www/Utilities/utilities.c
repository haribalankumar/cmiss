/*
 *  Utilities.c
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fopen(), fclose(), fgets(), NULL */

#include <stdlib.h>
    /* For: malloc(), free(), exit(), etc */

#include <string.h>
    /* For: strlen(), strtok(), etc */

#include <ctype.h>
    /* For: tolower() */

#include "htstring.h"
    /* For: htmlEscapedStrcpy(), htmlEscapedStrlen(), etc */

#include "utilities.h"
    /* For: out public interface */



/* -- Small Utilities ---------------------------------------------------------
*/ 

void forceLower( char *string )
  {
  char *ptr;

  ptr = string;

  while( *ptr != '\0' )
    {
    *ptr = tolower( *ptr );
    ptr++;
    }
  }



/* -- HTML Output methods -----------------------------------------------------
*/

BOOLEAN htmlString( char *string )
  {
  char bodyPrefix[7];
  
  strncpy( bodyPrefix, string, 6 );
  bodyPrefix[6] = '\0'; /* make sure of termination */
  forceLower( bodyPrefix );

  if( strcmp( bodyPrefix, "<html>" ) == 0 )  
    return True;
  else
    return False;
  }


#define Normal 0
#define Newline 1

void outputString( FILE *output, char *string )
  {
  char *newString;
  int   state;
  char *ptr;
  char *endPtr;
  int   length;
  BOOLEAN paragraphPending;

  if( htmlString( string ) == True )
    {
    /* Output the bits between the HTML tags */
    ptr = string + 6;
    endPtr = strstr( "</HTML>", string );
    if( endPtr == NULL )
      endPtr = string + strlen( string );
    fwrite( ptr, 1, endPtr - ptr, output );
    }
  else
    {
    length = htmlEscapedStrlen( string );
    newString = malloc( length + 1 );
    htmlEscapedStrcpy( newString, string );

    paragraphPending = False;
    ptr = newString;
    state = Normal;
    while( *ptr != '\0' )
      {
       if( state == Normal )
        {
        if( *ptr == '\n' )
          state = Newline;
        else
          if( paragraphPending == True )
            {
            fprintf( output, "<P>" );
            paragraphPending = False;
            }
        }
      else
        {
        if( *ptr != '\t' || *ptr != ' ' ) state = Normal;
        if( *ptr == '\n' ) paragraphPending = True;
        }
     fputc( *ptr++, output );
      }
    }
  }

