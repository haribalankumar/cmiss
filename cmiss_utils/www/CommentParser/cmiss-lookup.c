/*
 *  CMISS Lookup.c
 *
 *    Some test code fr playing around with the database functions
 */


/* -- Include Directives ------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stderr */

#include <stdlib.h>
    /* For: malloc(), free(), NULL */

#include <string.h>
    /* For: strcpy() */

#include <ctype.h>
    /* For: isprint(), tolower() */

#include <sys/types.h>
#include <sys/stat.h>
    /* For: stat() */

#include "wwwpaths.h"
    /* For CMISS_WWW_ROOT LOOKUP_URLPATH EXAMPLES_URLPATH LOOKUP_PROGRAM */


#include "groups.h"
    /* For: GROUP, GROUPFILE, setGroupFilePos(), etc, etc */

#include "rules.h"
    /* For: RULES, newRules(), etc */

#include "rule-parse.h"
    /* For: parseRules() */

#include "cgi-decode.h"
  /* For: CGI handling functions */

#include "htstring.h"
    /* For: htmlEscapedStrcpy(), htmlEscapedStrlen() */

#include "utilities.h"
    /* For: outputString() */

#include "reducers.h"
    /* For: reduceVariable(), reduceModule(), reduceSubroutine(), etc */

#include "simple-parser.h"
    /* For: PARSESTATE, newParseState(), extract(), skip(), etc */

#include "regexp.h"
    /* For: regexp, regcomp(), regexec() */

#include "outputcommands.h"
    /* For  outputCommand(), outputCommand2() etc */

/* -- Module Datatypes --------------------------------------------------------
*/

typedef void GROUPOUTPUT_FUNCTION( FILE *output, GROUP *group );




/* -- Module Constants --------------------------------------------------------
*/

char *VariableBase   = CMISS_WWW_ROOT LOOKUP_URLPATH "/variables";
char *ModuleBase     = CMISS_WWW_ROOT LOOKUP_URLPATH "/modules";
char *SubroutineBase = CMISS_WWW_ROOT LOOKUP_URLPATH "/subroutines";
char *CommentBase    = CMISS_WWW_ROOT LOOKUP_URLPATH "/comments";
char *CommandBase    = CMISS_WWW_ROOT LOOKUP_URLPATH "/commands";
char *FunctionBase   = CMISS_WWW_ROOT LOOKUP_URLPATH "/functions";
char *BlockDataBase   = CMISS_WWW_ROOT LOOKUP_URLPATH "/blockdata";

char *CommonFilePath = CMISS_ROOT "/cm/source/";



/* -- Private Module Methods --------------------------------------------------
*/



/* -- Group Printing Routines -------------------------------------------------
*/

void outputVariable( FILE *output, GROUP *group )
  {
  TAG        *tag;

  if( group == NULL )
    return;
    
  fprintf( output, "<HR>" );

  fprintf( output, "<H2>%s</H2>", group->name );

  tag = lookupTag( group, "TYPE" );
  if( tag != NULL )
    {
    fprintf( output, "<B>%s</B>", tag->body );
    }

  tag = lookupTag( group, "SET-UP" );
  if( tag != NULL )
    {
    fprintf( output, " - <I>Set up in " );
    outputNameList( output, tag->body, CommaSeperatedList );
    fprintf( output, "</I>" );
    }
  else
    {
    /* fprintf( output, " - <I>Defined in file %s</I>", group->file ); */
    }

  tag = lookupTag( group, "DESCRIPTION" );
  if( tag != NULL )
    {
    fprintf( output, "<P>" );
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


void outputSubroutine( FILE *output, GROUP *group )
  {
  TAG        *tag;
  char       *routineName;
  PARSESTATE *parse;

  if( group == NULL )
    return;

  fprintf( output, "<HR>" );

  fprintf( output, "<H2>%s</H2>", group->name );

  tag = lookupTag( group, "DEFINED-IN" );
  if( tag != NULL )
    {
    fprintf( output, "<I>(Defined in Module " );
    fprintf( output, "<A HREF=\"%s?name=%s\">%s</A>", LOOKUP_PROGRAM, 
      tag->body, tag->body );
    fprintf( output, ")</I><P>" );
    }
    
  tag = lookupTag( group, "DESCRIPTION" );
  if( tag != NULL )
    {
    fprintf( output, "<P>" );
    outputString( output, tag->body );
    }

  tag = lookupTag( group, "CALLS" );
  if( tag != NULL )
    {
    fprintf( output, "<P><DL><DT><B>Calls:</B><DD>" );
    outputNameList( output, tag->body, CommaSeperatedList );
    fprintf( output, "</DL>" );
    }

  tag = lookupTag( group, "CALLED-FROM" );
  if( tag != NULL )
    {
    fprintf( output, "<P><DL><B>Called From:</B><DD>" );
    outputNameList( output, tag->body, CommaSeperatedList );
    fprintf( output, "</DL>" );
    }

  tag = lookupTag( group, "SEE-ALSO" );
  if( tag != NULL )
    {
    fprintf( output, "<P><B>See also:</B> " );
    outputNameList( output, tag->body, CommaSeperatedList );
    }

  tag = lookupTag( group, "DEFINED-IN" );
  if( tag != NULL )
    {
    parse = newParseState( group->name );
    if( parse == NULL )
      {
      /* Error reporting? */
      return;
      }

    skip( parse, " \t\n\r" );
    routineName = extract( parse, " \t\n\r" );
    
    fprintf( output, "<FORM METHOD=\"GET\""
      "ACTION=\""
      "routine-browser.cgi\">"
      "<INPUT NAME=\"routine\" TYPE=\"hidden\" VALUE=\"%s\">"
      "<INPUT TYPE=\"submit\" VALUE=\"View Routine\">"
      "</FORM>", routineName );
    
    free( routineName );
    disposeParseState( parse );
    }

  fprintf( output, "<P>" );
  }


void outputFunction( FILE *output, GROUP *group )
  {
  TAG        *tag;
  PARSESTATE *parse;
  char       *routineName;

  if( group == NULL )
    return;

  fprintf( output, "<HR>" );

  fprintf( output, "<H2>%s</H2>", group->name );

  tag = lookupTag( group, "TYPE" );
  if( tag != NULL )
    {
    fprintf( output, "<B>%s</B>", tag->body );
    }
  
  tag = lookupTag( group, "DEFINED-IN" );
  if( tag != NULL )
    {
    fprintf( output, "<I>(Defined in Module " );
    fprintf( output, "<A HREF=\"%s?name=%s\">%s</A>", LOOKUP_PROGRAM, 
      tag->body, tag->body );
    fprintf( output, ")</I><P>" );
    }
    
  tag = lookupTag( group, "DESCRIPTION" );
  if( tag != NULL )
    {
    fprintf( output, "<P>" );
    outputString( output, tag->body );
    }

  tag = lookupTag( group, "CALLS" );
  if( tag != NULL )
    {
    fprintf( output, "<P><DL><DT><B>Calls:</B><DD>" );
    outputNameList( output, tag->body, CommaSeperatedList );
    fprintf( output, "</DL>" );
    }

  tag = lookupTag( group, "CALLED-FROM" );
  if( tag != NULL )
    {
    fprintf( output, "<P><DL><B>Called From:</B><DD>" );
    outputNameList( output, tag->body, CommaSeperatedList );
    fprintf( output, "</DL>" );
    }

  tag = lookupTag( group, "SEE-ALSO" );
  if( tag != NULL )
    {
    fprintf( output, "<P><B>See also:</B> " );
    outputNameList( output, tag->body, CommaSeperatedList );
    }

  tag = lookupTag( group, "DEFINED-IN" );
  if( tag != NULL )
    {
    parse = newParseState( group->name );
    if( parse == NULL )
      {
      /* Error reporting? */
      return;
      }

    skip( parse, " \t\n\r" );
    routineName = extract( parse, " \t\n\r" );
  
    fprintf( output, "<FORM METHOD=\"GET\" ACTION=\"routine-browser.cgi\">"
	    "<INPUT NAME=\"routine\" TYPE=\"hidden\" VALUE=\"%s\">"
	    "<INPUT TYPE=\"submit\" VALUE=\"View Routine\">"
	    "</FORM>", routineName );

    free( routineName );
    disposeParseState( parse );
    }

  fprintf( output, "<P>" );
  }


void outputBlockData( FILE *output, GROUP *group )
  {
  TAG        *tag;
  PARSESTATE *parse;
  char       *routineName;

  if( group == NULL )
    return;

  fprintf( output, "<HR>" );

  fprintf( output, "<H2>%s</H2>", group->name );

  tag = lookupTag( group, "TYPE" );
  if( tag != NULL )
    {
    fprintf( output, "<B>%s</B>", tag->body );
    }
  
  tag = lookupTag( group, "DESCRIPTION" );
  if( tag != NULL )
    {
    fprintf( output, "<P>" );
    outputString( output, tag->body );
    }

  tag = lookupTag( group, "SEE-ALSO" );
  if( tag != NULL )
    {
    fprintf( output, "<P><B>See also:</B> " );
    outputNameList( output, tag->body, CommaSeperatedList );
    }

  tag = lookupTag( group, "DEFINED-IN" );
  if( tag != NULL )
    {
    fprintf( output, "<I>(Defined in Module " );
    fprintf( output, "<A HREF=\"%s?name=%s\">%s</A>", LOOKUP_PROGRAM, 
      tag->body, tag->body );
    fprintf( output, ")</I><P>" );
    }

  parse = newParseState( group->name );
  if( parse == NULL )
    {
    /* Error reporting? */
    return;
    }

  skip( parse, " \t\n\r" );
  routineName = extract( parse, " \t\n\r" );

  fprintf( output, "<FORM METHOD=\"GET\""
    "ACTION=\""
    "routine-browser.cgi\">"
    "<INPUT NAME=\"routine\" TYPE=\"hidden\" VALUE=\"%s\">"
    "<INPUT TYPE=\"submit\" VALUE=\"View Routine\">"
    "</FORM>", routineName );

  free( routineName );
  disposeParseState( parse );

  fprintf( output, "<P>" );
  }



void outputComment( FILE *output, GROUP *group )
  {
  TAG *tag;

  if( group == NULL )
    return;

  fprintf( output, "<HR>" );

  fprintf( output, "<H2>%s</H2>", group->name );

  tag = lookupTag( group, "DESCRIPTION" );
  if( tag != NULL )
    {
    fprintf( output, "<P>" );
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


void outputModule( FILE *output, GROUP *group )
  {
  TAG *tag;
  int  numberOfTags;
  int  i;
  PARSESTATE *state;
  char *name;
  char *description;
#if 0
  char *moduleName;
  char *moduleNametmp;
#endif
 
  if( group == NULL )
    return;

  fprintf( output, "<HR>" );

  fprintf( output, "<H2>%s</H2>", group->name );

  tag = lookupTag( group, "DESCRIPTION" );
  if( tag != NULL )
    {
    outputString( output, tag->body );
    fprintf( output, "<P>" );
    }

  tag = lookupTag( group, "SEE-ALSO" );
  if( tag != NULL )
    {
    fprintf( output, "<P><B>See also:</B> " );
    outputNameList( output, tag->body, CommaSeperatedList );
    }

#if 0 /* There is not usually a file associated with a module */
  state = newParseState( group->name );
  if( state == NULL )
    {
    /* Error reporting? */
    return;
    }
  skip( state, " \t\n\r" );
  moduleName = extract( state, " \t\n\r" );
  moduleNametmp = malloc( strlen( moduleName ) +5 );
  reduceModule( moduleNametmp, moduleName);
  strcat(moduleNametmp,".f");
  fprintf( output, "<FORM METHOD=\"GET\""
    "ACTION=\""
    "routine-browser.cgi\">"
    "<INPUT NAME=\"module\" TYPE=\"hidden\" VALUE=\"%s\">"
    "<INPUT TYPE=\"submit\" VALUE=\"View Source\">"
    "</FORM>", moduleNametmp );
  free( moduleName );
  free( moduleNametmp );
  disposeParseState( state );
#endif

  tag = lookupTag( group, "ROUTINES" );
  if( tag != NULL )
    {
      fprintf( output, "<P><B>Routines:</B> " );
      outputNameList( output, tag->body, UnorderedBulletList );
    }

  numberOfTags = getNumberOfTags( group, "ROUTINE" );
  if( numberOfTags != 0 )
    {
    for( i = 1; i <= numberOfTags; i++ )
      {
      tag = lookupNthTag( group, "ROUTINE", i );
      if( tag != NULL )
        {
        state = newParseState( tag->body );
        if( state == NULL )
          {
          /* out of memory, leave */
          return;
          }
        skip( state, " \t\n\r" );
        name = extract( state, " \t\n\r" );
        if( name != NULL )
          {
          fprintf( output, "<A HREF=\"%s?name=%s\">%s</A> - ", LOOKUP_PROGRAM, name, name );
          free( name );
          }
        skip( state, " \t\n\r" );
        description = extract( state, "" );
        if( description != NULL )
          {
          outputString( output, description );
          free( description );
          }
        fprintf( output, "<BR>" );
        
        disposeParseState( state );
        }
      }
    }

  fprintf( output, "<P>" );
  }



/* -- Module Methods ----------------------------------------------------------
*/
#if 0
static void lookupKey( const char *key, const char *baseName,  
  REDUCTION_FUNCTION *reduce, GROUPOUTPUT_FUNCTION *output, int *count )
  {
  char       word[128];
  ENTRY      entry;
  FILEPOS    position;
  DATABASE  *database;
  RULES     *rules;
  GROUPFILE *file;
  GROUP     *group;
  char      *fileName;

  fileName = malloc( strlen( baseName ) + 7 );
  if( fileName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in outputGroup()\n" );
    return;
    }

  sprintf( fileName, "%s", baseName );
  database = openDatabase( fileName, ReadOnlyDatabase );
  if( database == NULL )
    {
    /* Some sort of error reported here? */
    return;
    }

  reduce( word, key );

  entry.key = word;
  entry.keyLength = strlen( word ) + 1;

  readEntry( database, &entry );

  if( entry.content == NULL )
    {
    return;
    }
    
  memcpy( &position, entry.content, entry.contentLength );
      /* -- Why is this call nessessary, and why wont a cast work here?
            Answers on the back of an envelope please */
 
  sprintf( fileName, "%s.list", baseName );
  file = newGroupFile( fileName );
  if( file == NULL )
    {
    /* Error reporting here? Error condition passed back to our caller? */
    return;
    }

  rules = newRules();
  if( rules == NULL )
    {
    /* Error reporting here? Error condition passed back to our caller? */
    return;
    }

  sprintf( fileName, "%s.rules", baseName );
  parseRuleFile( rules, fileName );

  setGroupFilePosition( file, position );

  group = getGroup( file, rules );

  if( group != NULL )
    {
    output( stdout, group );
    (*count)++;
    }
  }
#endif /* 0 */

void mungeRegularExpression( char *destination, const char *source )
  {
  const char *s;
  char       *d;
  const char *current;
  const char *prev;

  /* What we want to do is force all input to be lower case, bracket
     the who expression with "^" and "$" to make sure we only get full
     matches, and replace all "*"'s with ".*"'s where the "*" isn't
     preceeded by a close bracket (for now - other RE forms before a
     "*" may follow) */ 

  /* Ok - well the above worked quite well.
     Now for plan B:
       o  Keep the lowercasing code
       o  Keep the inplicit "^" and "$"
       o  Put a "." in front of every "?+*" *unless* the previous char
          was one of "]." (and maybe ")"?)
       o  Remove whitespace? (needed for comments? commands?)
  */

  /* Well that works too. Will wonders never cease? OK - I'll leave it
     at that until I get people complaining about it not being the same
     as the emacs search facility. */

  s = source;
  d = destination;
  prev = " ";  /* Anything for now */

  *d++ = '^';

  while( *s != '\0' )
    {
    current = s;

    if( isupper( *s ) )
      {
      *d++ = tolower( *s ); s++;
      }
    else if( *s == '*' || *s == '+' || *s == '?' )
      {
      if( *prev == ']' || *prev == ')' || *prev == '.' )
        {
        *d++ = *s++;
        }
      else
        {
        *d++ = '.';  /* Put an implicit "." before the wildcard */
        *d++ = *s++;
        }
      }
    else
      {
      if( isspace( *s ) )
        {
        s++;
        }
      else
        {
        *d++ = *s++;
        }
      }

    prev = current;
    }

  *d++ = '$';
  *d = '\0';
  }


void lookupRegularExpression( char *regularExpression, char *baseName,
  REDUCTION_FUNCTION *reduce, GROUPOUTPUT_FUNCTION *output, char *searchTag,
  int *count )
  {
  char      *fileName;
  char      *searchString;
  regexp    *compiledRegularExpression;
  GROUPFILE *file;
  RULES     *rules;
  GROUP     *group;
  TAG       *tag;
  char      *name;
  
  searchString = malloc( strlen( regularExpression ) * 2 + 3 );
  if( searchString == NULL )
    { 
    fprintf( stderr, "Memory allocation failure in "
      "lookupRegularExpression() [1]\n" );
    return;
    }

  mungeRegularExpression( searchString, regularExpression );

  compiledRegularExpression = regcomp( searchString );
  if( compiledRegularExpression == NULL )
    {
    fprintf( stderr, "regcomp() returned NULL in "
      "lookupRegularExpression()\n" );
    return;
    }

  fileName = malloc( strlen( baseName ) + 7 );
  if( fileName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in "
      "lookupRegularExpression() [1]\n" );
    return;
    }

  sprintf( fileName, "%s.list", baseName );

  file = newGroupFile( fileName );
  if( file == NULL )
    {
    fprintf( stderr, "newGroupFile() returned NULL in "
      "lookupRegularExpression()\n" );
    return;
    }

  rules = newRules();
  if( rules == NULL )
    {
    fprintf( stderr, "newRules() returned NULL in "
      "lookupRegularExpression()\n" );
    return;
    }

  sprintf( fileName, "%s.rules", baseName );
  parseRuleFile( rules, fileName );
      /* -- No error checking... */

  group = getGroup( file, rules );
  while( group != NULL )
    {
    tag = lookupTag( group, searchTag );
    if( tag == NULL )
      {
      fprintf( stderr, "lookupTag() returned NULL in "
        "lookupRegularExpression()\n" );
      return;
      }

    if( ! ( name = malloc( strlen( tag->body ) + 1 ) ) )
      {
	fprintf( stderr,
"Memory allocation failure in lookupRegularExpression at line %d in %s\n",
		 __LINE__, __FILE__ );
	return;
      }

    reduce( name, tag->body );

    if( regexec( compiledRegularExpression, name ) == 1 )
      {
      output( stdout, group );
      (*count)++;
      }
    
    free( name );
    discardGroup( file );
    group = getGroup( file, rules );
    }
  }


BOOLEAN nameContainsRegularExpression( char *name ) 
  {
  char *p;
  
  /* Scan the string looking for anything in the list ".*+([" */
  p = name;
  while( *p )
    {
    if( *p == '.' || *p == '*' || *p == '+' || *p == '(' || *p == '[' )
      return True;
    p++;
    }

  return False;
  }



/* -- Common Browser ---------------------------------------------------------- 
*/

void truncateLine( char *line )
  {
  int i;

  for( i = strlen( line ) - 1; i >= 0; i-- )
    {
    if( line[i] > 31 && line[i] != '\n' && line[i] != '\r' && line[i] != '\t' &&        line[i] != ' ' )
      {
      line[i+1] = '\0';
      break;
      }
    else
      {
      if( i == 0 )
        line[i] = '\0';
      }
    }
  }


void lookupCommonFile( char *name, int *count )
  {
  char *extension;
  char fullName[512];
  char moduleName[512];
  struct stat stat_buf;

  if( strchr( name, '/' ) != NULL )
    {
    /* "/" is illegal in the filenames... */
    return;
    }
/* KAT 3Nov99: Only searching for file if it has the appropriate extension */

  extension = strchr( name, '.' );
  if( extension == NULL )
    {
      return;
    }
  extension++;
  if( strcmp( extension, "cmn" ) != 0 && strcmp( extension, "inc" ) != 0
     && strcmp( extension, "f" ) != 0 && strcmp( extension, "c" ) != 0 )
    {
      return;
    }

  sprintf( moduleName, "%s", name );
  sprintf( fullName, "%s%s", CommonFilePath, name );

  /* This is a bit of a simplistic check really */
/*  if( strchr( name, '.' ) == NULL )
    {
    sprintf( moduleName, "%s.cmn", name );
    sprintf( fullName, "%s%s.cmn", CommonFilePath, name );
    }
  else
    {
    sprintf( moduleName, "%s", name );
    sprintf( fullName, "%s%s", CommonFilePath, name );
    }
*/

  if( 0 != stat( fullName, &stat_buf ) )
    {
    /* Can't stat the file. Perhaps the user had typed in something that wasn't
       a common file name. Punt back to our caller. */
    return;
    }
  
  printf( "<H2>%s</H2><P>", moduleName );

  printf("<FORM METHOD=\"GET\""
    "ACTION=\""
    "routine-browser.cgi\">"
    "<INPUT NAME=\"module\" TYPE=\"hidden\" VALUE=\"%s\">"
    "<INPUT TYPE=\"submit\" VALUE=\"View Source\">"
    "</FORM>", moduleName );

  printf( "</PRE>\n" );

  /* Indicate that we output something for our caller */
  (*count)++;
  }



/* -- Misc --------------------------------------------------------------------
*/


char *Form = 
  "<HR><P>"
  "Please type in a full name used within CMISS, or a regular expression "
  "that matches one or more names within CMISS (see <A HREF="
  "\""
  LOOKUP_URLPATH
  "/cmiss-browser-help.html\">"
  "help on regular expression searches</A>).\n"
  "<FORM METHOD=\"GET\" "
  "ACTION=\"" LOOKUP_PROGRAM "\">"
  "CMISS Name: <INPUT NAME=\"name\" SIZE=30> "
  "<INPUT TYPE=\"submit\" VALUE=\"CMISS Lookup\">"
  "</FORM>\n";

void printForm( void )
  {
  printf( Form );
  }



/* -- Test program entry point ------------------------------------------------
*/

#define Header \
  "Content-type: text/html\n\n" \
  "<HTML><HEAD><TITLE>CMISS Viewer</TITLE></HEAD><BODY>" \
  "<H1>CMISS Viewer</H1>" \
  "\n"

#define Footer \
  "<HR>" \
  "<ADDRESS>" \
  "<A HREF=\"" \
  CMISS_WWW_URLPATH \
  "/help/index_programmer.html\">CMISS Programmer Help</A>" \
  "</ADDRESS></BODY></HTML>\n"
 
int main( int argc, char *argv[] )
  {
  CGI  *cgi;
  char *name;
  int  count;

  printf( Header ); 
  fflush( stdout );  /* Hack to make sure stderrr doesn't get seen first */

  cgi = getCGIEnvironment( argc, argv );
  if( cgi == NULL )
    {
    printf( "<HR><P><B>ERROR:</B> getCGIEnvironment() failed in main()</P>"
            "<P>Email e.poor@auckland.ac.nz <I>immediately</I></P>\n" );
    }
  else
    {
    printForm();

    name = lookupString( cgi, "name" );

    count = 0;

    if( name != NULL )
      {
/*  KAT 10Sep99: Doing a complete search every time in case a key occurs more
 than once in the database */
/*      lookupKey( name, VariableBase, reduceVariable, outputVariable,
        &count );
      lookupKey( name, ModuleBase, reduceModule, outputModule, &count );
      lookupKey( name, SubroutineBase, reduceSubroutine, outputSubroutine,
        &count );
      lookupKey( name, CommentBase, reduceComment, outputComment, &count );
      lookupKey( name, CommandBase, reduceCommand, outputCommand, &count );
      lookupKey( name, FunctionBase, reduceFunction, outputFunction, &count );
      lookupKey( name, BlockDataBase, reduceBlockData, outputBlockData, &count );

      if( count == 0 )  / * Nothing matched exactly * /
        {
        / * OK - attempt 2: try assuming that the name is really a common
           file * /


         lookupCommonFile( name, &count );

        if( count == 0 )
          {
          / * Still nothing - do a full RE search * /
      if( nameContainsRegularExpression( name ) == True )
	{ */
      lookupRegularExpression( name, VariableBase, reduceVariable,
              outputVariable, "VARIABLE", &count );
      lookupRegularExpression( name, ModuleBase, reduceModule,
              outputModule, "MODULE", &count );
      lookupRegularExpression( name, SubroutineBase, reduceSubroutine,
              outputSubroutine, "SUBROUTINE", &count );
      lookupRegularExpression( name, CommentBase, reduceComment,
              outputComment, "COMMENT", &count );
      lookupRegularExpression( name, CommandBase, reduceCommand,
              outputCommand, "COMMAND", &count );
      lookupRegularExpression( name, FunctionBase, reduceFunction,
              outputFunction, "FUNCTION", &count );
      lookupRegularExpression( name, BlockDataBase, reduceBlockData,
              outputBlockData, "BLOCKDATA", &count );

      lookupCommonFile( name, &count );
/*          }   
          }
        } */

      if( count == 0 )  /* No search matched at all */
        {
        printf( "<HR><H2>Search Failed</H2>The name \"%s\" could not be "
          "found in the CMISS databases. Please check your spelling and "
          "try again.<P>", name );
        }
      }
    }

  printf( Footer );

  return 0;
  }
 
