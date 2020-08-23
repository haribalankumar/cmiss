/*
 * Edit File.c
 *
 *    Brings up a file for the user to edit.
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
    /* For: JOB, getJob(), printJobHeading(), etc */

#include "directory.h"
    /* For: FILEINFO, DIRECTORY, newDirectory() etc */

#include "webcmiss.h"


/* -- HTML Output -------------------------------------------------------------
*/

void printEditForm( FILE *output, char *fullFileName, char *leafFileName,
  USER *user )
  {
  FILESTATE *file;
  char      *line;

  file = newFileState( fullFileName );
  if( file == NULL )
    {
    fprintf( output, "Error" );
    return;
    }

  fprintf( output, "<CENTER>" );

  fprintf( output, "<TABLE BORDER=\"0\" WIDTH=\"80%%\"><TR>" );

  fprintf( output, "<TD ALIGN=\"CENTER\" COLSPAN=\"3\">" );

  printf( "<FORM METHOD=\"POST\" ACTION=\"" SAVE_FILE "\">\n" );
  printf( "<INPUT NAME=\"file\" TYPE=\"hidden\" VALUE=\"%s\">", leafFileName );
  printf( "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">", user->name );
  printf( "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">", 
    user->password );

  printf( "<TEXTAREA ROWS=\"32\" COLS=\"68\" NAME=\"text\">\n" );
  line = getLine( file );
  while( line != NULL )
    {
    printf( "%s", line );
    discardLine( file );
    line = getLine( file );
    }
  printf( "</TEXTAREA>" );

  fprintf( output, "</TD>" );

  fprintf( output, "<TR>" );
  fprintf( output, "<TD ALIGN=\"LEFT\">" );
  fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"Save File\">" );
  fprintf( output, "</TD>" );

  fprintf( output, "<TD ALIGN=\"CENTER\">" );
  fprintf( output, "<INPUT TYPE=\"RESET\" VALUE=\"Undo all changes\">" );
  fprintf( output, "</FORM>" );
  fprintf( output, "</TD>" ); 

  fprintf( output, "<TD ALIGN=\"RIGHT\">" );
  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" LIST_PROBLEM "\">" );
  fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">",
    user->name );
  fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">", 
    user->password );
  fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"Cancel\">" );
  fprintf( output, "</FORM>" );
  fprintf( output, "</TD>" );
  
  fprintf( output, "</TABLE>" );

  fprintf( output, "</CENTER>" );

  disposeFileState( file );
  }



/* -- Module Methods ----------------------------------------------------------
*/

void editFile( CGI *cgi, FILE *output )
  {
  USER      *user;
  JOB       *job;
  char      *fullFileName;
  char      *fileName;

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
    fprintf( stderr, "getJob() returned NULL in editFile()\n" );
    disposeUser( user );
    return;
    }

  fileName = lookupString( cgi, "file" );
  if( fileName == NULL )
    {    
    fprintf( stderr, "No file to edit in editFile()\n" );
    return;
    } 

  if( strchr( fileName, '/' ) != NULL || strstr( "..", fileName ) != NULL )
    {
    fprintf( stderr, "Illegal character in filename in editFile()\n" );
    return; 
    }
  
  fullFileName = malloc( strlen( job->inputDirectoryName ) +
    strlen( fileName ) + 2 ); /* + '/' + '\0' */
  if( fullFileName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in editFile()\n" );
    return;
    }

  sprintf( fullFileName, "%s/%s", job->inputDirectoryName, fileName );

  fprintf( output, "<HR NOSHADE><P>" );
  printJobHeading( output, user, job );
  fprintf( output, "<P>" );
 
  printEditForm( output, fullFileName, fileName, user );
  }



/* -- Program Entry Point -----------------------------------------------------
*/

#define Header \
  "Content-type: text/html\n\n" \
  "<HTML><HEAD><TITLE>File Editor</TITLE></HEAD>" \
  "<BODY BGCOLOR=\"#ffffff\">" \
  "<H1>File Editor</H1>" \
  "\n"

#define Footer \
  "<HR NOSHADE>\n"

int main( int argc, char *argv[] )
  {
  CGI       *cgi;

  printf( Header ); 
  fflush( stdout );  /* Hack to make sure stderr doesn't get seen first */

  cgi = getCGIEnvironment( argc, argv );
  if( cgi == NULL )
    {
    printf( "<HR><P><B>ERROR:</B> getCGIEnvironment() failed in main()</P>"
            "<P>Email "WEBMASTER" <I>immediately</I></P>\n" );
    return 0;
    }

  editFile( cgi, stdout );

  printf( Footer );

  return 0;
  }

