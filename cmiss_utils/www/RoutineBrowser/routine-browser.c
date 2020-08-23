/*
 *  View Routine.c
 *
 *    CGI-Bin program to load and view a portion of CMISS source
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
/* For: printf() */

#include <errno.h>

#include <stdlib.h>
/* For: exit(), getenv() */

#include <string.h>
/* For: strcmp(), strtok(), etc */

#include <ctype.h>
/* For: toupper() */

#include "wwwpaths.h"
    /* For CMISS_WWW_ROOT, CMISS_WWW_URLPATH */

#include "htstring.h"
/* For: htmlEscapedStrlen(), htmlEscapedStrcpy() */

#include "simple-parser.h"
/* For: PARSESTATE, newParseState(), skip(), prefix() etc */

#include "../CGIDecode/cgi-decode.h"
/* For: CGI handling functions */

#include <sys/types.h>
#include <unistd.h>
/* added bt M. Stevenson for getpid */

/* -- Datatypes ---------------------------------------------------------------
*/

typedef struct
{
  char *fileName;
  char *moduleName;
  char *routineName;
  long  startByte;
  long  startLine;
  long  endByte;
  long  endLine;
}
ENTRY;



/* -- Test setup --------------------------------------------------------------
*/

void setupTest( void )
{
  extern int putenv(char *);

  putenv( "SERVER_SOFTWARE=NCSA/1.3" );
  putenv( "SERVER_NAME=www.esc.auckland.ac.nz" );
  putenv( "GATEWAY_INTERFACE=CGI/1.1" );
  putenv( "SERVER_PROTOCOL=HTTP/1.0" );
  putenv( "SERVER_PORT=80" );
  putenv( "REQUEST_METHOD=GET" );
  putenv( "SCRIPT_NAME=/cgi-bin/edouards-test" );
  putenv( "QUERY_STRING=routine=FEM" );
  putenv( "REMOTE_HOST=esu6.aukuni.ac.nz" );
  putenv( "REMOTE_ADDR=130.216.5.106" );
}



/* -- Module Constants --------------------------------------------------------
*/

char *CommonFilePath = CMISS_ROOT "/cm/source/";
char *allow_ip_prefix1 = "130.216.208.";
char *allow_ip_prefix2 = "130.216.218.";


/* -- Private Methods ---------------------------------------------------------
*/

void forceUpper( const char *source, char *destination )
{
  const char *s = source;
  char       *d = destination;

  while( *s )
  {
    *d = toupper( *s );
    d++;
    s++;
  }
}



/*
 *  This is a mess!
 */
ENTRY *lookupEntry( char *entryFile, char *targetModule, char *targetRoutine )
{
  FILE  *file;
  char   line[256];
  char  *name;
  char  *module;
  char  *routine;
  ENTRY *entry;

  if( targetModule != NULL )
  {
    if( *targetModule == '\0' )
    {
      targetModule = NULL;
    }
    else
    {
      forceUpper( targetModule, targetModule ); 
    }
  }

  forceUpper( targetRoutine, targetRoutine );

  file = fopen( entryFile, "rt" );
  if( file == NULL )
    { 
      fprintf( stderr, "Failed to open file %s: %s\n",
	       entryFile, strerror(errno) );
      return NULL;
    }
  
  fgets( line, 256, file );
  while( !feof( file ) )
  {
    name = strtok( line, " " );
    module = strtok( NULL, " " );
    routine = strtok( NULL, " " );

    if( module != NULL && routine != NULL )
    {
      forceUpper( routine, routine );
      forceUpper( module, module );

      if( strcmp( routine, targetRoutine ) == 0 &&
        ( targetModule == NULL || strcmp( module, targetModule ) == 0 ) )
      {
        entry = malloc( sizeof( ENTRY ) );
       
        entry->fileName = malloc( strlen( name ) + 1 );
        entry->moduleName = malloc( strlen( module ) + 1 );
        entry->routineName = malloc( strlen( routine ) + 1 );
        /* -- WARNING: No Error Checking Done */

        strcpy( entry->fileName, name );
        strcpy( entry->moduleName, module );
        strcpy( entry->routineName, routine );

        sscanf( strtok( NULL, " " ), " %ld ", &(entry->startByte) );
        sscanf( strtok( NULL, " " ), " %ld ", &(entry->startLine) );
        sscanf( strtok( NULL, " " ), " %ld ", &(entry->endByte) );
        sscanf( strtok( NULL, " " ), " %ld ", &(entry->endLine) );
  
        return entry;
      }
    }
    
    fgets( line, 256, file );
  }

  /* Not found... */
  return NULL;
}

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

#define EntryFile CMISS_WWW_ROOT CM_URLPATH "/routines/routine-index"

BOOLEAN subroutineMatch( char *routineName, char *line )
{
  PARSESTATE *state;
  BOOLEAN     foundFlag;

  state = newParseState( line );
  if( state == NULL )
  {
    fprintf( stderr, "newParseState() returned NULL in subroutineMatch()\n" );
    return False;
  }

  skip( state, " \t\r\n" );

  if( empty( state ) == True )
  {
    foundFlag = False;
  }
  else if( prefix( state, "SUBROUTINE" ) == True )
  {
    skipWord( state, "SUBROUTINE" );
    skip( state, " \t\r\n" );
    if( prefix( state, routineName ) == True )
    {
      foundFlag = True;
    }
    else
    {
      foundFlag = False;
    }
  }
  else if( prefix( state, "BLOCK DATA" ) == True )
  {
    skipWord( state, "BLOCK DATA" );
    skip( state, " \t\r\n" );
    if( prefix( state, routineName ) == True )
    {
      foundFlag = True;
    }
    else
    {
      foundFlag = False;
    }
  }
  else
  {
    skipUntil( state, " \t\r\n" );
    skip( state, " \t\r\n" );
    if( prefix( state, "FUNCTION" ) == True )
    {
      skipWord( state, "FUNCTION" );
      skip( state, " \t\r\n" );
      if( prefix( state, routineName ) == True )
      {
        foundFlag = True;
      }
      else
      {
        foundFlag = False;
      }
    }
    else
    {
      foundFlag = False;
    }
  }

  disposeParseState( state );

  return foundFlag;
}


void viewRoutine( char *moduleName, char *routineName )
{
  size_t      length;
  ENTRY      *theEntry;
  FILE       *theFile;
  char       *theBuffer;
  char       *escapedBuffer;
  char        line[256];
  char        upperCaseLine[256];
  char        escapedline[256];
  long        linesRead;
  BOOLEAN     endFound;
  PARSESTATE *state;
  
  forceUpper( routineName, routineName );
  /* -- Actually this is rather bad practice - the place where
  routineName was originally defined should be in charge of
  making sure it's upper case. Doing it here on this side of
  function call means that it isn't const, and that
  the caller may well get back a string that differed from
  that passed in... */

  theEntry = lookupEntry( EntryFile, moduleName, routineName );
  if( theEntry != NULL )
  {
    printf( "<H1>%s</H1>\n", theEntry->routineName );

    theFile = fopen( theEntry->fileName, "rt" );
    if( theFile == NULL )
      { 
	printf( "<p><b>ERROR:</b> Failed to open file %s: %s</p>\n",
		theEntry->fileName, strerror(errno) );
	return;
      }
    fseek( theFile, theEntry->startByte, SEEK_SET );
 
    /* read a single line to see if the index for this routine is still
    valid */
    fgets( line, 256, theFile );
    fseek( theFile, theEntry->startByte, SEEK_SET );
    forceUpper( line, line );
    if( subroutineMatch( routineName, line ) == True )
    {
      length = theEntry->endByte - theEntry->startByte;
      theBuffer = malloc( length + 1 );
    
      fread( theBuffer, 1, length, theFile );

      theBuffer[length] = '\0';

      length = htmlEscapedStrlen( theBuffer );
      escapedBuffer = malloc( length + 1 );
      htmlEscapedStrcpy( escapedBuffer, theBuffer );

      printf( "<P><I>(see lines %ld to %ld in module %s)</I></P>\n", 
        theEntry->startLine, theEntry->endLine, theEntry->moduleName );

#ifndef NOT_SETGID
      /* Perl colorisation script won't run */
      printf( "<PRE>\n" );
      puts( escapedBuffer );
      printf( "</PRE>\n" );

#else /* not setgid: can use Perl colorisation script */
      {
	/* Addition by M. Stevenson to add syntax hilighting to code */
	FILE *tempfile;
	char filename[30]; /* more than enought to hold /tmp/pid.tmp */
	char perlscript[80];
	/* this should be dynamic but 2048 should be bug enough */
	/* need the size because html'fied code can be big */
	char buff[2048]; 
	pid_t pid;

      /* added by M.S name of temp file */
      pid = getpid();
      sprintf(filename,"/tmp/color%d.tmp",(int) pid);

      /* open temp file */
      if( !(tempfile = fopen(filename,"w")) ) printf("\nCould not open %s\n",filename);
      /* write sub to temp file */
      fprintf(tempfile,"%s",escapedBuffer);
      fclose(tempfile);

      /* run perl colorisation script on file keep name same */
      sprintf(perlscript,"/www/Groups/Bioengineering/CMISS/scripts/color_source.pl %s",filename);
      system(perlscript);
      /* reopen temp file now with colorised html code */
      if( !(tempfile = fopen(filename,"r")) ) printf("\nCould not open %s\n",filename);
      /* printf( "<P>\n" ); */
      while( fgets(buff, 2047, tempfile) != NULL ){
        puts( buff );
      }

      unlink(filename);
      }
      /* end of M.S changes */
#endif

      fclose( theFile );

      free( theBuffer );
      free( escapedBuffer );
    }
    else
    {
      if( theEntry->startByte < 60000 ) 
      fseek( theFile, 0, SEEK_SET );
      else
      fseek( theFile, theEntry->startByte - 60000, SEEK_SET );

      linesRead = 0;
      fgets( line, 256, theFile );
      forceUpper( line, upperCaseLine );
      linesRead++;
      while( !feof( theFile ) &&
        linesRead < (theEntry->endLine - theEntry->startLine + 2000) )
      {
        if( subroutineMatch( routineName, upperCaseLine ) == True )
        {
          /* OK. We've found it. */
          printf( "<PRE>" );

          htmlEscapedStrcpy( escapedline, line );
          printf( "%s", escapedline );
          
          fgets( line, 256, theFile );
          endFound = False;
          while( endFound != True )
          {
            state = newParseState( line );
            if( prefix( state, " " ) == True )
            {
              shiftChars( state, 6 );
              if( prefix( state, "END" ) == True )
              {
                skipWord( state, "END" );
                skip( state, " \t\n\r" );
                if( empty( state ) == True || prefix( state, "!" ) == True )
                endFound = True;
              }
            }
            disposeParseState( state );
            htmlEscapedStrcpy( escapedline, line );
            printf( "%s", escapedline );

            fgets( line, 256, theFile );
            if( feof( theFile ) )
            endFound = True;
          }
          printf( "</PRE>" );
          break;
        }

        fgets( line, 256, theFile );
        if( !feof( theFile ) )
        {
          forceUpper( line, upperCaseLine );    
          linesRead++;
        }
      }
    }
  }
  else
  {
    printf( "<P>The previous search for the routine <B>%s</B> failed"
      " (does the routine exist?)</P>\n", routineName );
  }
}

#define MaxLineLength 256

void viewModule( char *moduleName )
{
  FILE *file;
  char *temp;
  char line[MaxLineLength];
  char fullName[512];


  sprintf( fullName, "%s%s", CommonFilePath, moduleName );
  file = fopen( fullName, "rt" );
  if( file == NULL )
  {
    /* Can't open it. Perhaps the user had typed in something that wasn't
    a common file name. Punt back to our caller. */
    return;
  }

  printf( "<P><B>Filename:</B> <I>%s</I></P>", moduleName );

  printf( "<PRE>\n" );

  fgets( line, MaxLineLength, file );
  while( !feof( file ) )
  {
    truncateLine( line );
    temp = malloc( htmlEscapedStrlen( line ) + 1 );
    htmlEscapedStrcpy( temp, line );
    printf( "%s\n", temp );
    free( temp );
    fgets( line, MaxLineLength, file );
  }

  printf( "</PRE>\n" );

}

/* -- The Program Entry Point -------------------------------------------------
*/

#define Header \
"Content-type: text/html\n\n" \
"<HTML>" \
"<HEAD>" \
"<TITLE>CMISS Routine Browser</TITLE>" \
"</HEAD>" \
"<BODY>" \
"<H1>CMISS Routine Browser</H1>\n"

#define OldForm \
"<HR>" \
"<FORM METHOD=\"GET\" ACTION=\"/Groups/Bioengineering/CMISS/scripts/routine-browser.cgi\">" \
"<P><B>Please input the CMISS routine you wish to view:</B></P>" \
"Routine Name: <INPUT NAME=\"routine\" SIZE=30><BR>" \
"<B>Then press either:</B> <INPUT TYPE=\"submit\" VALUE=\"View Routine\"> or " \
"<INPUT TYPE=\"reset\" VALUE=\"Clear Fields\">" \
"</FORM>" \
"<HR>"
/*"Module Name: <INPUT NAME=\"module\" SIZE=30> <I>(optional)</I><BR><BR>" \ */

#define Form \
"<HR><P>" \
"Please type in a full name used within CMISS, or a regular expression " \
"that matches one or more names within CMISS (see <A HREF=" \
  "\"" \
  LOOKUP_URLPATH \
  "/cmiss-browser-help.html\">" \
"help on regular expression searches</A>).\n" \
"<FORM METHOD=\"GET\" " \
"ACTION=\"" LOOKUP_PROGRAM_URLPATH "\">" \
"CMISS Name: <INPUT NAME=\"name\" SIZE=30> " \
"<INPUT TYPE=\"submit\" VALUE=\"CMISS Lookup\">" \
"</FORM>" \
"<HR>\n"


#define Footer \
"<ADDRESS>" \
"<A HREF=\"" \
CMISS_WWW_URLPATH \
"/help/index_programmer.html\">CMISS Programmer Help</A> / " \
"Routine Browser" \
"</ADDRESS></BODY></HTML>\n"


int main( int argc, char *argv[] )
{
  CGI   *cgi;
  char  *moduleName;
  char  *routineName;
  char  *remote_addr;

#if defined( TEST )
  setupTest();
#endif

  printf( Header );

  printf( Form );

  cgi = getCGIEnvironment( argc, argv );

  remote_addr = getenv( "REMOTE_ADDR" );

  if( cgi == NULL )
  {
    printf( "<HR><P><B>ERROR:</B> getCGIEnvironment() failed in main()</P>"
      "<P>Email cmisshelp@esu1.auckland.ac.nz</P>\n" );
  }
  else if( remote_addr == NULL )
    {
      printf( "<p><b>ERROR:</b> REMOTE_ADDR unknown.</p>\n" );
    }
  else if( 0 != strncmp( remote_addr, allow_ip_prefix1,
			 strlen(allow_ip_prefix1) )
	   && 0 != strncmp( remote_addr, allow_ip_prefix2,
			    strlen(allow_ip_prefix2) )
	   )
    {
      printf( "<p>Source code is not public.\n</p>" );
    }
  else
  {
    moduleName = lookupString( cgi, "module" );
    routineName = lookupString( cgi, "routine" );

    if( routineName != NULL && *routineName != '\0' )
    {
      viewRoutine( moduleName, routineName );
      printf( "<HR>" );
    }
    if( moduleName != NULL && *moduleName != '\0' )
    {
      viewModule(moduleName);
      printf( "<HR>" );
    }
  }
  printf( Footer );

  disposeCGI( cgi );

  return 0;
}
