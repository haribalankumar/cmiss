/*
 *  View File.c
 *
 *    Views a file in the Output directory.
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
  /* For: fprintf(), stderr, etc */

#include <stdlib.h>
  /* For: malloc(), free(), NULL */

#include <string.h>
  /* For: strlen(), strchr(), strstr(), strcpy() */

#include "cgi-decode.h"
  /* For: CGI, getCGIEnvironment(), lookupString(), etc */

#include "user-managment.h"
    /* For: USER, USERLIST, readUserList(), passwordCheck(), etc */

#include "get-user.h"
    /* For: getUser() */

#include "job.h"
    /* For: JOB, getJob(), printJobHeading(), etc */

#include "directory.h"
    /* For: FILEINFO, DIRECTORY, newDirectory() etc */

#include "webcmiss.h"

/* -- Module Constants --------------------------------------------------------
*/


/* -- Private Module Methods --------------------------------------------------
*/

void displayFile( FILE *output, char *fullFileName )
  {
  FILE *file;
  char  line[512];

  file = fopen( fullFileName, "rt" );
  if( file == NULL )
    {
    fprintf( stderr, "fopen() failed in displayFile()\n" );
    return;
    }

  fprintf( output, "<TABLE WIDTH=\"90%%\" BGCOLOR=\"#eeeeee\" BORDER=\"0\">\n" );
  fprintf( output, "<TR><TD>\n" );
  fprintf( output, "<PRE>\n" );
 
  fgets( line, 512, file );
  while( !feof( file ) )
    {
    fprintf( output, "%s", line );
    fgets( line, 512, file );
    }

  fprintf( output, "</PRE>\n" );
  fprintf( output, "</TD></TR></TABLE>\n" );
  }


void printReturnButton( FILE *output, USER *user, char *return_script, char *directory )
  {
  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"%s\">\n", return_script );
  fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->name );
  fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->password );
  fprintf( output, "<INPUT NAME=\"directory\" TYPE=\"hidden\" VALUE=\"%s\">\n", directory );
  fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"Return\">\n" );
  fprintf( output, "</FORM>\n" );
  }


void viewFile( CGI *cgi, FILE *output )
  {
  USER      *user;
  JOB       *job;
  char      *fullFileName;
  char      *dirName;
  char      *fileName;
	char      *directory;
  char      *return_script;

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

  if(!(directory = lookupString(cgi,"directory")))
    {
    directory = "input";
    }
  if (strcmp(directory,"root") == 0) {
    dirName = "";
  }
  else if (strcmp(directory,"input") == 0) {
    dirName = "Input";
  }
  else {
    dirName = "Output";
  }

  if(!(return_script = lookupString(cgi,"return")))
    {
    return_script = LIST_PROBLEM;
    }

  if( strchr( fileName, '/' ) != NULL || strstr( "..", fileName ) != NULL )
    {
    fprintf( stderr, "Illegal character in filename in editFile()\n" );
    return; 
    }
  
  fullFileName = malloc( strlen( user->homeDirectory ) + strlen( dirName ) +
    strlen( fileName ) + 3 ); /* + '/' + '\0' */
  if( fullFileName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in editFile()\n" );
    return;
    }

  sprintf( fullFileName, "%s/%s/%s", user->homeDirectory, dirName, fileName );

  fprintf( output, "<HR NOSHADE><P>\n" );
  printJobHeading( output, user, job );
  fprintf( output, "<P>\n" );

  fprintf( output, "<CENTER>\n" );
  fprintf( output, "<TABLE BORDER=\"0\" WIDTH=\"65%%\">\n" );
  fprintf( output, "<TR><TD ALIGN=\"LEFT\">\n" );
  fprintf( output, "<FONT SIZE=\"+1\">You are viewing the file "
    "<B>%s</B></FONT>\n", fileName );
  fprintf( output, "</TD><TD ALIGN=\"RIGHT\">\n" );
  printReturnButton( output, user, return_script, directory );
  fprintf( output, "</TD></TR></TABLE>\n" );
  fprintf( output, "</CENTER>\n" );
 
  fprintf( output, "<CENTER>\n" );
  displayFile( output, fullFileName );
  fprintf( output, "</CENTER>\n" );
  }



/* -- Program Entry Point -----------------------------------------------------
*/

#define Header \
  "Content-type: text/html\n\n" \
  "<HTML><HEAD><TITLE>File Viewer</TITLE></HEAD>" \
  "<BODY BGCOLOR=\"#ffffff\">" \
  "<H1>File Viewer</H1>" \
  "\n"

#define Footer \
  "<HR NOSHADE>\n"

int main( int argc, char *argv[] )
  {
  CGI *cgi;

  printf( Header ); 
  fflush( stdout );  /* Hack to make sure stderr doesn't get seen first */

  cgi = getCGIEnvironment( argc, argv );
  if( cgi == NULL )
    {
    printf( "<HR><P><B>ERROR:</B> getCGIEnvironment() failed in main()</P>"
            "<P>Email "WEBMASTER" <I>immediately</I></P>\n" );
    return 0;
    }

  viewFile( cgi, stdout );

  return 0;
  }
