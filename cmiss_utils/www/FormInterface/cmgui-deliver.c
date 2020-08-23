/*
 *  CMGUI Deliver.c
 *
 *    Returns exnode/exelem files, plus a .com command file
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "get-user.h"
    /* For: getUser() */

#include "job.h"
    /* For: JOB, getJob(), etc */

#include "cgi-decode.h"
    /* For: CGI handling functions */

#include "directory.h"
    /* For: DIRECTORY, FILEINFO, etc */

#include "webcmiss.h"


/* -- Module Constants --------------------------------------------------------
*/




/* -- Private Module Methods --------------------------------------------------
*/


void deliverFile( CGI *cgi )
  {
  USER *user;
  JOB  *job;
  char *directory;
  char *fileName;
  char *command;
  char *tarcom = "tar cf - *.ex* ";

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
    fprintf( stderr, "getJob() returned NULL in deliverFile()\n" );
    return;
    }
  
  directory = lookupString( cgi, "directory" );
  if( directory == NULL )
    {    
    directory = "output";
    }

  fileName = lookupString( cgi, "file" );
  if( fileName == NULL )
    {
    fprintf( stderr, "No file name specified in deliverFile()\n" );
    return;
    }
  
  if( strchr( fileName, '/' ) != NULL || strstr( "..", fileName ) != NULL )
    {
    fprintf( stderr, "Illegal character in filename in deliverFile()\n" );
    return; 
    }

  if (chdir(job->outputDirectoryName))
    {
    fprintf( stderr, "Error in chdir in deliverFile()\n" );
    return; 
    }

  command = malloc( strlen( fileName ) + strlen( tarcom ) + 2 );
  if( command == NULL )
    {
    fprintf( stderr, "Memory allocation failure in deliverFile()\n" );
    return;
    }

  sprintf( command, "%s %s", tarcom, fileName );

  system( command );

  free(command);

  forwardUserDir( LIST_PROBLEM, user, directory );
  }



/* -- Program Entry Point -----------------------------------------------------
*/

#define Header \
  "Content-type: application/x-cmgui\n\n"

int main( int argc, char *argv[] )
  {
  CGI       *cgi;

  printf( Header ); 
  fflush( stdout );  /* Hack to make sure stderr doesn't get seen first */

  cgi = getCGIEnvironment( argc, argv );
  if( cgi == NULL )
    {
    return 0;
    }

  deliverFile( cgi );

  return 0;
  }
