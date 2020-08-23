/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stderr */

#include <stdlib.h>
    /* For: NULL, malloc(), free() */

#include <string.h>
    /* For: strlen() */

#include <ctype.h>
    /* For: isprint(), tolower() */

#include <sys/types.h>
#include <sys/stat.h>
    /* For: stat() */

#include "wwwpaths.h"
    /* For CMISS_WWW_ROOT LOOKUP_URLPATH EXAMPLES_URLPATH LOOKUP_PROGRAM */





#include "groups.h"
    /* For: GROUP, GROUPLIST, newGroupList(), etc */

#include "rules.h"
    /* For: RULES, newRules(), etc */

#include "rule-parse.h"
    /* For: parseRules() */

#include "utilities.h"
    /* For: outputString() */


#include "simple-parser.h"
    /* For: PARSESTATE, newParseState(), extract(), skip(), etc */


/* -- Module Datatypes --------------------------------------------------------
*/

typedef enum
  {
  CommaSeperatedList,
  UnorderedBulletList
  }
LISTFORMAT;

/* -- HTML Output Helper Routines ---------------------------------------------
*/

/*
 *  The following code assumes that the name passed in will be in the
 *  form "2a3", which means that the resulting example file will
 *  be found in ".../example_2/example_2a/example_2a3/example_2a3.com"
 *
 *  This is a bit messy, but never mind.
 */
char *generateExampleDirectoryName( char *examples_path, char *exampleName )
  {
  char *fullName;
  char *index;
  int   length;
  int   fullLength;
  int   i;

  /* Handle backwards compatibility (welcome to the rest of your life) */
  index = strstr( exampleName, "example_" );
  if( index != NULL )
    {
    exampleName = index + strlen( "example_" ); 

    while( *index != '\0' && *index != '.' ) 
      index++;

    *index = '\0'; 
    }

  length = strlen( exampleName );
  fullLength = 9 * ( length + 1 ) + length * length; /* Over-allocates */
    /* Actually the previous line over-allocated before, but now that
       I've commented out the strcat line below, it vastly over-allocates.
       This is in need of a fix */

  fullName = malloc( strlen( examples_path ) + 1 + fullLength );

  strcpy( fullName, examples_path );

  for( i = 0; i < length; i++ )
    {
    /*strcat( fullName, "example_" );*/
    strncat( fullName, exampleName, i+1 );
    strcat( fullName, "/" );
    }

  /*
  strcat( fullName, "example_" );
  strcat( fullName, exampleName );
  strcat( fullName, ".com" );
  */

  return fullName;
  }


void outputNameList( FILE *output, char *string, LISTFORMAT format )
  {
  PARSESTATE *parse;
  char       *name;
  char       *urlEncodedName;
  BOOLEAN     firstName;

  parse = newParseState( string );
  if( parse == NULL )
    {
    /* Error reporting here? */
    return; 
    }

  firstName = True;

  while( empty( parse ) != True )
    {
    skip( parse, " \t,\n\r" );
    if( empty( parse ) != True )
      {
      name = extract( parse, ",\n\r" );
      
      switch( format )
        {
        case CommaSeperatedList:
          if( firstName == True )
            firstName = False;
          else
            fprintf( output, ", " );
          break;
  
        case UnorderedBulletList:
          if( firstName == True )
            {
            fprintf( output, "<UL>" );
            firstName = False;
            }
          fprintf( output, "<LI>" );
          break;
        }

      urlEncodedName = malloc( urlEncodedStrlen( name ) + 1 );
      if( urlEncodedName == NULL )
        {
        /* Report Error? */
        free( name );
        break;
        }

      urlEncodedStrcpy( urlEncodedName, name );

      fprintf( output, "<A HREF=\"%s?name=%s\">%s</A>", LOOKUP_PROGRAM,
        urlEncodedName, name );

      free( urlEncodedName );
      free( name );
      }

    skip( parse, " \t,\n\r" );
    }

  switch( format )
    {
    case CommaSeperatedList:
      break;

    case UnorderedBulletList:
      fprintf( output, "</UL>" );
      break;
    }


  disposeParseState( parse );
  }


void outputExampleList( FILE *output, char *string, LISTFORMAT format )
  {
  PARSESTATE *parse;
  char       *name;
  char       *path;
  char       *urlEncodedName;
  BOOLEAN     firstName;

  parse = newParseState( string );
  if( parse == NULL )
    {
    /* Error reporting here? */
    return; 
    }

  firstName = True;

  while( empty( parse ) != True )
    {
    skip( parse, " \t,\n\r" );
    if( empty( parse ) != True )
      {
      name = extract( parse, ",\n\r" );
      
      switch( format )
        {
        case CommaSeperatedList:
          if( firstName == True )
            firstName = False;
          else
            fprintf( output, ", " );
          break;
  
        case UnorderedBulletList:
          if( firstName == True )
            {
            fprintf( output, "<UL>" );
            firstName = False;
            }
          fprintf( output, "<LI>" );
          break;
        }
      
      path = generateExampleDirectoryName( EXAMPLES_URLPATH, name );
      if( path == NULL )
      {   
        fprintf( stderr, "generateExampleDirectoryName() returned NULL in "
          "changeExample()\n" );
        return;
      }

      urlEncodedName = malloc( urlEncodedStrlen( name ) + 1 );
      if( urlEncodedName == NULL )
        {
        /* Report Error? */
        free( name );
        break;
        }

      urlEncodedStrcpy( urlEncodedName, name );
      fprintf( output, "<A HREF=\"%sgenerated/test_page.html\">%s</A>", path,
        urlEncodedName);

      free( urlEncodedName );
      free( name );
      }

    skip( parse, " \t,\n\r" );
    }

  switch( format )
    {
    case CommaSeperatedList:
      break;

    case UnorderedBulletList:
      fprintf( output, "</UL>" );
      break;
    }


  disposeParseState( parse );
  }


void outputCommand( FILE *output, GROUP *group )
  {
  TAG        *tag;
  int         numberOfTags;
  int         i;
  PARSESTATE *state;
  char       *string;

  if( group == NULL )
    return;

  fprintf( output, "<HR>" );

  /*
   *  There was a request to remove this heading...
   *
  fprintf( output, "<H2>" );
  outputString( output, group->name );
  fprintf( output, "</H2>" );
   *
   */

  fprintf( output, "<DL>" );
  fprintf( output, "<DT>" );
  fprintf( output, "<B>" );
  outputString( output, group->name );
  fprintf( output, "</B>" );

  numberOfTags = getNumberOfTags( group, "PARAMETER" );
  for( i = 1; i <= numberOfTags; i++ )
    {
    tag = lookupNthTag( group, "PARAMETER", i );
    if( tag != NULL )
      {
      fprintf( output, "<DD>" );

      /* Ok - If the line is more than one line, then everything on
         line two and on is a description of the paramter. */
      state = newParseState( tag->body );
      if( state != NULL )
        {
        /* this is the parameter line */
        string = extract( state, "\n\r" );
        fprintf( output, "<B>" );
        outputString( output, string );
        fprintf( output, "</B>" );
        free( string );

        /* this is the optional description */
        skip( state, "\r\n" );
        if( empty( state ) != True )
          {
          string = extract( state, "" );  /* -- "" means end of string */
          fprintf( output, "<DL><DD>" );
          outputString( output, string );
          fprintf( output, "</DL><P>" );
          free( string );
          }

        disposeParseState( state );
        }
      else
        {
        outputString( output, tag->body );
        }
      }
    }
  fprintf( output, "</DL></B>" );

  tag = lookupTag( group, "DESCRIPTION" );
  if( tag != NULL )
    {
    fprintf( output, "<P>" );
    outputString( output, tag->body );
    }

  tag = lookupTag( group, "SEE-EXAMPLE" );
  if( tag != NULL )
    {
    fprintf( output, "<P><B>See example:</B> " );
    outputString( output, tag->body );
    }

  tag = lookupTag( group, "SEE-ALSO" );
  if( tag != NULL )
    {
    fprintf( output, "<P><B>See also:</B> " );
    outputNameList( output, tag->body, CommaSeperatedList );
    }

  fprintf( output, "<P>" );
  }

void outputCommand2( FILE *output, GROUP *group )
{
  TAG        *tag;
  int         numberOfTags;
  int         i;
  PARSESTATE *state;
  char       *string;

  if( group == NULL )
  return;



  fprintf( output, "<DT>" );

  numberOfTags = getNumberOfTags( group, "PARAMETER" );
  for( i = 1; i <= numberOfTags; i++ )
  {
    tag = lookupNthTag( group, "PARAMETER", i );
    if( tag != NULL )
    {
      fprintf( output, "<DD>" );

      /* Ok - If the line is more than one line, then everything on
      line two and on is a description of the paramter. */
      state = newParseState( tag->body );
      if( state != NULL )
      {
        /* this is the parameter line */
        string = extract( state, "\n\r" );
         outputString( output, string );
        free( string );


        disposeParseState( state );
      }
      else
      {
        outputString( output, tag->body );
      }
    }
  }
  fprintf( output, "</DL></B>" );

  tag = lookupTag( group, "DESCRIPTION" );
  if( tag != NULL )
    {
      fprintf( output, "<P>" );
      outputString( output, tag->body );
    }
      fprintf( output, "<P>" );


}

void outputCommand3( FILE *output, GROUP *group )
  {
  TAG        *tag;
  int         numberOfTags;
  int         i;
  PARSESTATE *state;
  char       *string;

  if( group == NULL )
    return;

  fprintf( output, "<DL>" );
  fprintf( output, "<DT>" );
  fprintf( output, "<B>" );
  outputString( output, group->name );
  fprintf( output, "</B>" );

  tag = lookupTag( group, "DESCRIPTION" );
  if( tag != NULL )
    {
      fprintf( output, "<DD>" );
      outputString( output, tag->body );
    }

  numberOfTags = getNumberOfTags( group, "PARAMETER" );
  for( i = 1; i <= numberOfTags; i++ )
    {
    tag = lookupNthTag( group, "PARAMETER", i );
    if( tag != NULL )
      {
      fprintf( output, "<DD>" );

      /* Ok - If the line is more than one line, then everything on
         line two and on is a description of the paramter. */
      state = newParseState( tag->body );
      if( state != NULL )
        {
        /* this is the parameter line */
        string = extract( state, "\n\r" );
        fprintf( output, "<B>" );
        outputString( output, string );
        fprintf( output, "</B>" );
        free( string );

        /* this is the optional description */
        skip( state, "\r\n" );
        if( empty( state ) != True )
          {
          string = extract( state, "" );  /* -- "" means end of string */
          fprintf( output, "<DL><DD>" );
          outputString( output, string );
          fprintf( output, "</DL><P>" );
          free( string );
          }

        disposeParseState( state );
        }
      else
        {
         outputString( output, tag->body );
        }
      }
    }
  fprintf( output, "</DL></B>" );

/*   tag = lookupTag( group, "DESCRIPTION" ); */
/*   if( tag != NULL ) */
/*     { */
/*     fprintf( output, "<P>" ); */
/*     outputString( output, tag->body ); */
/*     } */

  tag = lookupTag( group, "SEE-ALSO" );
  if( tag != NULL )
    {
    fprintf( output, "<P>See also: " );
    outputNameList( output, tag->body, CommaSeperatedList );
    }

  tag = lookupTag( group, "SEE-EXAMPLE" );
  if( tag != NULL )
    {
    fprintf( output, "<P>See example: " );
    outputExampleList( output, tag->body, CommaSeperatedList );
    }

  fprintf( output, "<P>" );
  }
