/*
 *  Parameter Job.c
 *  
 *    Job control for the web HPC parameter system: start and stop jobs.
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), fflush(), stdout, stderr, etc */

#include <stdlib.h>
    /* For: malloc(), free(), NULL */

#include <string.h>
    /* For: strlen(), strcpy() */

#include "user-managment.h"
    /* For: USER */

#include "get-user.h"
    /* For: getUser() */

#include "job.h"
    /* For: JOB */

#include "host.h"
    /* For: HOST, HOSTITEM etc */

#include "border.h"
    /* For: getBorderHeader etc */

#include "cgi-decode.h"
    /* For: CGI handling functions */

/* testing, testing, 1, 2, 3 ... */
#include <pwd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <unistd.h>
#include <signal.h>
#include <errno.h>
#include <unistd.h>
#include <netinet/in.h>
#include <limits.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <fcntl.h>
#include "transfer.h"
    /* For: Socket communication */


/* -- Module Constants --------------------------------------------------------
*/
#include "webcmiss.h"
#define BUFFER_SIZE     1024
#define CGI_BOUNDARY    "cmisscgiboundary"

/* -- Job delivery code ----------------------------------------------------------
*/

void deliverFile( CGI *cgi, FILE *output, char *example )
  {
    USER *user;
    JOB  *job;
    char *fileName;
    char *command;
    char *tarcom = "tar cf - *.ex* ";
    char *buff;
    int   len;

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
  
    fileName  = "param_view.com";
  
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

    fprintf( output, "Content-type: application/x-cmgui\n\n" ); 
    fflush( output );

    sprintf( command, "%s %s", tarcom, fileName );
    system( command );
    free( command );

    /* Forward onward to the example page */
    len  = strlen( example ) + 27;
    buff = malloc( len );
    if ( buff == NULL )
    {
      writeError( "malloc fails in paramJob()\n" );
      return;
    }

    sprintf( buff, "directory=input&example=%s", example );
    forwardUserTag( PARAM_EXAMPLE, user, buff );
    free( buff );

    return;
  }



/* -- Job submission code --------------------------------------------------------
*/

void printOKButton( FILE *output, USER *user, char *example )
  {
    fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" PARAM_EXAMPLE "\">" );
    fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">", user->name );
    fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">", user->password );
    fprintf( output, "<INPUT NAME=\"example\" TYPE=\"hidden\" VALUE=\"%s\">", example );
    fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"OK\">" );
    fprintf( output, "</FORM>" );
  }


void submitJobCGI( CGI *cgi, FILE *output, char *threadstring, char *hostnamestring, 
		   char *versionname, HOSTLIST *hostList, char *example )
  {
    USER *user;
    JOB  *job;
    HOST *host;
    int   i,len;
    char  param_name[16];
    char *buff;
    char *value_string;


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
      paramBorderWriteError( output, user, "getJob() returned NULL in submitJobCGI()\n" );
      disposeUser( user );
      return;
    }

    host = lookupHost( hostnamestring, hostList );
    if( host == NULL )
    {
      paramBorderWriteError( output, user, "getHost() returned NULL in submitJobCGI()\n" );
      disposeHost( host );
      return;
    }

    if (chdir(job->inputDirectoryName))
    {
      paramBorderWriteError( output, user, "submitJobCGI: cannot cd to input directory\n" );
      return;
    }

    /* Update the input files: extract the parameter values */
    i = 1;
    while( i )
    {
      sprintf( param_name, "param_value_%d", i );
      value_string = lookupString( cgi, param_name );

      if ( value_string != NULL ) {
        buff = malloc( 16 + strlen( value_string ) );
        if ( buff == NULL ) {
          writeError( "malloc fails in paramJob()\n" );
          return;
        }

        sprintf( buff, "PARAM_VALUE_%d=%s", i, value_string );
        free( value_string );
        putenv( buff );

        i++;
      }
      else {
        i = 0;
      }
    }

    /* Run the script on the file */
    system( "/bin/sh ./param.cgi" );

    /* Run the job */
    paramBorderGenerateHeader( output, user, False, NULL, NULL );
    fprintf( output, "<CENTER>\n");
    switch ( submitJob( output, versionname, threadstring, host, job, user) )
    {
      case 0:
      {
        fprintf( output, "Your job has been submitted to %s.<BR> You will be sent email "
          "informing you that your job has finished its run.", hostnamestring );
        fflush( output );
        break;
      }
      case 1:
      {
        fprintf( output, "<B>Error:</B> Your previously submitted CMISS job has "
          "yet to complete on the University of Auckland's High Performance "
          "Computer. Please try again later. (Remember that you will be sent "
          "email when a job has been successfully completed.)" );
        fflush( output );
        break;
      }
      default:
      {
        fprintf( output, "Error in the job submission. " );
        fflush( output );
        break;
      }
    }
    fprintf( output, "</CENTER>\n");
    fprintf( output, "<P>\n");

    fprintf( output, "<CENTER>" );
    printOKButton( output, user, example );
    fprintf( output, "</CENTER>\n" );

    borderGenerateFooter( output );
    fflush( output );

    /* Forward onward to the example page */
    len  = strlen( example ) + 27;
    buff = malloc( len );
    if ( buff == NULL )
    {
      writeError( "malloc fails in paramJob()\n" );
      return;
    }

    sprintf( buff, "directory=input&example=%s", example );
    forwardUserTag( PARAM_EXAMPLE, user, buff );
    free( buff );

    return;
  }


/* -- Job abort code --------------------------------------------------
*/

void abortJobCGI( CGI *cgi, FILE *output, HOSTLIST *hostList, char *example )
  {
    USER *user;
    JOB  *job;
    HOST *host;
    HOSTITEM *hostItem;
    int   len;
    char *buff;

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
      paramBorderWriteError( output, user, "getJob() returned NULL in paramAbortJob()\n" );
      disposeUser( user );
      return;
    }

    host = lookupHost( job->host, hostList );
    if( host == NULL )
    /* Clear and abort all hosts */
    {
      hostItem = hostList->firstHostItem;
      while (hostItem != NULL)
      {
	abortJob( output, hostItem->host, job, user, 1);
	hostItem = hostItem->nextHostItem;
      }
    }
    else
    /* Abort on a named host */
    {
      /* rsh to hpc and do your business */
      if ( abortJob( output, host, job, user, 1) )
      {
	paramBorderWriteError( output, user, "abort() returned NULL in paramAbortJob()\n" );
	return;
      }
    }

    /* Forward onward to the example page */
    len  = strlen( example ) + 27;
    buff = malloc( len );
    if ( buff == NULL )
    {
      paramBorderWriteError( output, user, "malloc fails in paramJob()\n" );
      return;
    }

    sprintf( buff, "directory=input&example=%s", example );
    forwardUserTag( PARAM_EXAMPLE, user, buff );
    free( buff );
    
    return;
  }


/* -- Job view code --------------------------------------------------
*/

void viewJobCGI( CGI *cgi, FILE *output, char *example )
  {
    USER *user;
    JOB  *job;
    int   len;
    char *buff;

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
      paramBorderWriteError( output, user, "getJob() returned NULL in paramViewJob()\n" );
      disposeUser( user );
      return;
    }

    /* Forward onward to the example page */
    len  = strlen( example ) + 27;
    buff = malloc( len );
    if ( buff == NULL )
    {
      paramBorderWriteError( output, user, "malloc fails in paramJob()\n" );
      return;
    }

    sprintf( buff, "directory=input&example=%s", example );
    forwardUserTag( PARAM_EXAMPLE, user, buff );
    free( buff );
    
    return;
  }


/* -- Program Entry Point -----------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
    CGI *cgi;
    USER *user;
    char *threadsstring;
    char *hostnamestring;
    char *versionname;
    char *stagestring;
    char *example;
    int   stage;
    HOSTLIST *hostList;

    if ( ( versionname = malloc(BUFFER_SIZE) ) == NULL)
    {
      printf( "<HR><P><B>ERROR:</B> malloc() failed in main()</P>\n");
      return 0;
    }

    cgi = getCGIEnvironment( argc, argv );
    if( cgi == NULL )
    {
      printf( "<HR><P><B>ERROR:</B> getCGIEnvironment() failed in main()</P>"
	      "<P>Email "WEBMASTER" <I>immediately</I></P>\n" );
      return 0;
    }

    user = getUser( cgi );
    hostList = newHostList();
    readHostList( hostList, ESU_MASTER_HOST_LIST );

    /* What problem are we working on? */
    example = lookupString( cgi, "example" );
    if( example == NULL )
    {
      paramBorderWriteError( stdout, user, "No example name passed to paramJob()\n" );
      return 0;
    }

    if(!(stagestring = lookupString(cgi,"stage")))
    {
      stage = 0;
    }
    else
    {
      stage = atoi(stagestring);
    }

    /* If anything has been canceled, then set to the default */
    switch (stage)
    {
      /* Stage 0 -- submit job */
      case 0:
      {
	hostnamestring = "hpc1";
	threadsstring  = "1";
	versionname    = "cmo_32";
	submitJobCGI( cgi, stdout, threadsstring, hostnamestring, versionname, hostList, example );
	break;
      }
	
      /* Stage 10 -- abort job */
      case 10:
      {
	abortJobCGI( cgi, stdout, hostList, example );
	break;
      }

      /* Stage 20 -- view job */
      case 20:
      {
	deliverFile( cgi, stdout, example );
	/* viewJobCGI( cgi, stdout, example ); */
	break;
      }

      default:
      {
	paramBorderWriteError( stdout, user, "Invalid stage number" );
	break;
      }
    } 
    
    return 0;
  }
