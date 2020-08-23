/*
 *  Save File.c
 *
 *    Saves a file (i.e. FORM TEXTAREA data) back to disk
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: printf() */

#include <stdlib.h>
    /* For: NULL, malloc(), free() */

#include <string.h>
    /* For: strlen() */

#include "cgi-decode.h"
    /* For: CGI handling functions */

#include "user-managment.h"
    /* For: USER, USERLIST, readUserList(), passwordCheck(), etc */

#include "get-user.h"
    /* For: getUser() */

#include "job.h"
    /* For: JOB, getJob(), etc */

#include "directory.h"
    /* For: DIRECTORY, FILEINFO, etc */

#include "webcmiss.h"


/* -- Module Constants --------------------------------------------------------
*/ 




/* -- Module Data Structures --------------------------------------------------
*/




/* -- Method Prototypes -------------------------------------------------------
*/




/* -- Support Code ------------------------------------------------------------
*/ 

#ifdef UNUSED
/*
 *  dup()
 */
static char *newStringCopy( const char *source )
  {
  char *tempString;

  tempString = malloc( strlen( source ) + 1 );
  if( tempString == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newStringCopy()\n" );
    return NULL;
    }

  strcpy( tempString, source );
  
  return tempString;
  }
#endif


/*
 *  System Specific (i.e. this only work under UNIX) 
 */
char *newLeafString( const char *path )
  {
  char *tempString;
  char *ptr;

  /* Step 1. find leaf part of name */
  ptr = strrchr( path, '/'  );
  if( ptr == NULL )
    ptr = (char *) path;  /* path in const, ptr is not */
  else
    ptr += 1; /* skip the path element */

  /* Step 2. make a copy for our caller */
  tempString = malloc( strlen( ptr ) + 1 );
  if( tempString == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newLeafString()\n" );
    return NULL;
    }

  strcpy( tempString, ptr );
  
  return tempString;
  }






/* -- Save File ---------------------------------------------------------------
*/

void setHTMLMode( FILE *output )
  {
  fprintf( output, "Content-type: text/html\n\n" );
  fflush( output );  /* Hack to make sure stderr doesn't get seen first */
  }
 

void writeUnixTextFile( char *text, FILE *file )
  {
  char *p;
  int   last;
  /* Constants */
  const int Newline = 1;
  const int Alpha   = 2;

  p = text;
  last = Alpha;

  while( *p != '\0' )
    { 
    if( *p != '\n' && *p != '\r' )
      {
      fputc( *p, file );
      last = Alpha;
      }
    else
      {
      if( *p == '\n' )
        {
        fputc( '\n', file );
        last = Newline;
        if( *(p+1) == '\r' )
          p++;     /* Skip next char as well */
        }
      else if( *p == '\r' )
        {
        fputc( '\n', file );
        last = Newline;
        if( *(p+1) == '\n' )
          p++;    /* Skip next char as well */
        }
      }
    p++;
    }
  
  /* Hack to deal with daving a file with a missing final endline */
  if( last != Newline )
    fputc( '\n', file );  
  }


void saveFile( CGI *cgi, FILE *output )
  {
  USER      *user;
  JOB       *job;
  char      *fullFileName;
  char      *fileName;
  char      *text; 
  FILE      *file;

  user = getUser( cgi );
  if( user == NULL )
    {
    /* Not always an error, the user could simply have supplied an incorrect
       password. getUser() will report back those conditions to the user
       (one hopes...) */
    return;  
    }

  job = getJob( user );
  if( job == NULL )
    {
    fprintf( stderr, "getJob() returned NULL in revertFile()\n" );
    disposeUser( user );
    return;
    }

  fileName = lookupString( cgi, "file" );
  if( fileName == NULL )
    {    
    setHTMLMode( output );
    fprintf( stderr, "No filename to save to in saveFile()\n" );
    return;
    } 

  if( strchr( fileName, '/' ) != NULL || strstr( "..", fileName ) != NULL )
    {
    setHTMLMode( output );
    fprintf( stderr, "Illegal character in filename in saveFile()\n" );
    return; 
    }

  text = lookupString( cgi, "text" );
  if( text == NULL )
    {
    setHTMLMode( output );
    fprintf( stderr, "No file to save in saveFile()\n" );
    return;
    }

  fullFileName = malloc( strlen( job->inputDirectoryName ) +
    strlen( fileName ) + 2 );
  if( fullFileName == NULL )
    {
    setHTMLMode( output );
    fprintf( stderr, "Memory allocation failure in saveFile()\n" );
    return;
    }

  sprintf( fullFileName, "%s/%s", job->inputDirectoryName, fileName );

  /* remove file, then write a new one */
  if( remove( fullFileName ) != 0 )
    {
    setHTMLMode( output );
		printf( "fullFileName = %s\n", fullFileName );
    fprintf( stderr, "remove() failed in saveFile()\n" );
    free( fullFileName );
    return;
    }

  file = fopen( fullFileName, "wt" );
  if( file == NULL )
    {
    setHTMLMode( output );
		printf( "fullFileName = >%s<\n", fullFileName );
    fprintf( stderr, "fopen() returned NULL in saveFile()\n" );
    free( fullFileName );
    return;
    }

  writeUnixTextFile( text, file );

  fclose( file );

  /* If everything was successful */
	forwardUser( LIST_PROBLEM, user );
  }



/* -- Program Entry Point -----------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
  CGI       *cgi;

  cgi = getCGIEnvironment( argc, argv );
  if( cgi == NULL )
    {
    printf( "<HR><P><B>ERROR:</B> getCGIEnvironment() failed in main()</P>"
            "<P>Email "WEBMASTER" <I>immediately</I></P>\n" );
    return 0;
    }

  saveFile( cgi, stdout );

  return 0;
  }

