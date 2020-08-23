/*
 *  Job Control.c
 *  
 *    Job control for the web HPC server: start, stop and check status
 *    of jobs.
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

/* -- Job submission code --------------------------------------------------------
*/

void printOKButton( FILE *output, USER *user )
  {
    fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" JOB_CONTROL "\">" );
    fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">", user->name );
    fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">", user->password );
    fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"OK\">" );
    fprintf( output, "</FORM>" );
  }


void printJobAlreadyRunning( FILE *output, USER *user )
  {
    fprintf( output, "<FONT SIZE=\"+1\">" );
    fprintf( output, "<B>Error:</B> Your previously submitted CMISS job has "
	     "yet to complete on the University of Auckland's High Performance "
	     "Computer. Please try again later. (Remember that you will be sent "
	     "email when a job has been successfully completed.)" );
    fprintf( output, "</FONT>" );

    fprintf( output, "<CENTER>" );

    printOKButton( output, user );

    fprintf( output, "</CENTER>" );
  }


void submitJobCGI( FILE *output, CGI *cgi, char *threadstring, char *hostnamestring, 
		   char *versionname, HOSTLIST *hostList )
  {
    USER *user;
    JOB  *job;
    HOST *host;

    fprintf( output, "<H1>Job Submission</H1><P>" );
    fprintf( output, "<HR NOSHADE><P>" );

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
	fprintf( stderr, "getJob() returned NULL in submitJobCGI()\n" );
	disposeUser( user );
	return;
      }

    host = lookupHost( hostnamestring, hostList );
    if( host == NULL )
      {
	fprintf( stderr, "getHost() returned NULL in submitJobCGI()\n" );
	disposeHost( host );
	return;
      }

    if (chdir(job->inputDirectoryName))
      {
	fprintf( output, "submitJobCGI: cannot cd to \"%s\"\n",
		 job->inputDirectoryName);
	return;
      }

#ifndef TEST_CODE
    switch ( submitJob( output, versionname, threadstring, host, job, user) )
      {
      case 0:
	{
	  fprintf( output, "Your job has been submitted to %s. You will be sent email "
		   "informing you that your job has finished its run, and any "
		   "files created will be placed in the <B>Output</B> directory "
		   "accessible from the file listing on your CMISS web page.", hostnamestring );
	  fflush( output );
	  break;
	}
      case 1:
	{
	  printJobAlreadyRunning( output, user ); 
	  fflush( output );
	  return;
	}
      default:
	{
	  fprintf( output, "<HR NOSHADE><P>" );
	  fprintf( output, "Error in the job submission. " );
	  fflush( output );
	  break;
	}
      }
#else
    fprintf( output, "Your job has been submitted to %s using %s with %s threads\n",
	     hostnamestring,versionname,threadstring);
    fprintf( output, "The account %s is used\n",host->username);
    fprintf( output, "<P>" );
    fflush( output );
#endif
    fprintf( output, "<P>" );
  
    fprintf( output, "<CENTER>" );
    printOKButton( output, user );
    fprintf( output, "</CENTER>\n" );
  }

void writeHostForm( FILE *output, char *userstring, char *passwordstring, HOSTLIST *hostList )
  {
    HOSTITEM *hostItem;
    BOOLEAN first = True;

    fprintf( output, "<H1>CMISS Host Selection</H1><P>\n" );
    fprintf( output, "<HR NOSHADE><P>\n" );

    if ( hostList == NULL )
      {
	fprintf( output, "<p>No hosts defined. Cannot run any jobs.\n" );
	return;
      }

    fprintf( output, "<P><form method=\"POST\" action=\"" JOB_CONTROL "\">\n" );
    fprintf( output, "<input type=\"hidden\" name=\"user\" value=\"%s\">\n", userstring );
    fprintf( output, "<input type=\"hidden\" name=\"password\" value=\"%s\">\n", passwordstring );
    fprintf( output, "<input type=\"hidden\" name=\"stage\" value=\"2\">\n" );
    fprintf( output, "Select the host on which to run:<br>\n");
    fprintf( output, "<P>\n");
    fprintf( output, "<OL>\n");

    hostItem = hostList -> firstHostItem;
    while( hostItem != NULL )
      {
	fprintf( output, "<li><input type=\"radio\" name=\"hostname\" value=\"%s\"",
		 hostItem -> host -> hostname );
	if (first)
	  {
	    fprintf( output, " CHECKED" );
	    first = False;
	  }
	fprintf( output, ">%s (%s, %s cpu)\n", hostItem -> host -> hostname, 
		 hostItem -> host -> os, hostItem -> host -> ncpu );
	hostItem = hostItem -> nextHostItem;
      }
    fprintf( output, "</ol>\n" );

    fprintf( output, "<P>\n" );
    fprintf( output, "<ol>\n" );
    fprintf( output, "<input type=\"submit\" value=\"Submit\">\n" );
    fprintf( output, "<input type=\"reset\" value=\"Reset\">\n" );
    fprintf( output, "<input type=\"submit\" value=\"Cancel\" name=\"cancel\">\n" );
    fprintf( output, "</ol>\n" );
    fprintf( output, "</form></p>\n" );
}

void writeVersionForm( FILE *output, char *userstring, char *passwordstring,
		       char *hostnamestring, HOSTLIST *hostList )
  {
    int i,ncpu;
    HOST *host;

    fprintf( output, "<H1>CMISS Host Selection</H1><P>\n" );
    fprintf( output, "<HR NOSHADE><P>\n" );

    host = lookupHost( hostnamestring, hostList );
    ncpu = atoi(host->ncpu);

    fprintf( output, "<p><form method=\"POST\" action=\"" JOB_CONTROL "\">\n" );
    fprintf( output, "<input type=\"hidden\" name=\"user\" value=\"%s\">\n", userstring );
    fprintf( output, "<input type=\"hidden\" name=\"password\" value=\"%s\">\n", passwordstring );
    fprintf( output, "<input type=\"hidden\" name=\"hostname\" value=\"%s\">\n", hostnamestring );
    fprintf( output, "<input type=\"hidden\" name=\"stage\" value=\"3\">\n" );
    fprintf( output, "Select the version of CMISS to run on %s<br>\n", hostnamestring );
    fprintf( output, "<P>\n");

    fprintf( output, "Optimisation:\n");
    fprintf( output, "<ol>\n");
    {
      fprintf( output, "<li><input type=\"radio\" name=\"debug\" value=\"0\" CHECKED>Optimised\n");
      fprintf( output, "<li><input type=\"radio\" name=\"debug\" value=\"1\">Unoptimised\n");
    }
    fprintf( output, "</ol>\n");
    fprintf( output, "<P>\n");

    if ( strcmp( host->os, "irix_64" ) == 0 || strcmp( host->os, "aix" ) == 0 )
      {
	fprintf( output, "Pointer Size:\n");
	fprintf( output, "<ol>\n");
	{
	  fprintf( output, "<li><input type=\"radio\" name=\"pointer\" value=\"64\" CHECKED>64 bit\n");
	  fprintf( output, "<li><input type=\"radio\" name=\"pointer\" value=\"32\">32 bit\n");
	}
	fprintf( output, "</ol>\n");
	fprintf( output, "<P>\n");
      }
    else
      {
	fprintf( output, "<input type=\"hidden\" name=\"pointer\" value=\"32\">\n");			
      }

    if ( ncpu > 1 )
      {
	fprintf( output, "Serial/Parallel Executable:\n");
	fprintf( output, "<ol>\n");
	{
	  fprintf( output, "<li><input type=\"radio\" name=\"parallel\" value=\"0\" CHECKED>Single processor\n");
	  fprintf( output, "<li><input type=\"radio\" name=\"parallel\" value=\"1\">Parallel\n");
	}
	fprintf( output, "</ol>\n");
	fprintf( output, "<P>\n");

	fprintf( output, "Number of Processors:\n");
	{
	  fprintf( output, "<select name=\"threads\">\n");
	  fprintf( output, "<option selected>1\n");
	  for (i=2;i<=ncpu;i++)
	    {
	      fprintf( output, "<option>%d\n",i);
	    }
	  fprintf( output, "</select>\n");
	}
	fprintf( output, "<P>\n");
      }
    else
      {
	fprintf( output, "<input type=\"hidden\" name=\"parallel\" value=\"0\">\n");			
	fprintf( output, "<input type=\"hidden\" name=\"threads\" value=\"1\">\n");			
      }
  
    fprintf( output, "<ol>\n");
    fprintf( output, "<input type=\"submit\" value=\"Submit\">\n");
    fprintf( output, "<input type=\"reset\" value=\"Reset\">\n");
    fprintf( output, "<input type=\"submit\" value=\"Cancel\" name=\"cancel\">\n" );
    fprintf( output, "</ol>\n");
    fprintf( output, "</form></p>\n");
  }

int getVersion(CGI *cgi, char *executable, HOSTLIST *hostList)
  {
    int n;
    char *debug, *hostname, *parallel, *pointer;
    HOST *host;

    /* Hostname/OS */
    if(hostname = lookupString(cgi,"hostname"))
      {
	host = lookupHost( hostname, hostList );
	if ( host == NULL )
	  {
	    writeError("Cannot find host\n");
	    return -1;				
	  }
      }
    else
      {
	writeError("No host specified\n");
	return -1;
      }

    /* Optimisation level */
    if(debug = lookupString(cgi,"debug"))
      {
	if ( atoi(debug) )
	  {
	    strcpy(executable,"cm");
	    n = 2;
	  }
	else
	  {
	    strcpy(executable,"cmo");
	    n = 3;
	  }
      }
    else
      {
	writeError("No optimisation level selected\n");
	return -2;
      }

    /* Pointer size */
    if(pointer = lookupString(cgi,"pointer"))
      {
	if ( atoi(pointer) == 64 )
	  {
	    strcpy(executable+n,"_64");
	    n += 3;
	  }
	else
	  {
	    strcpy(executable+n,"_32");
	    n += 3;
	  }
      }
    else
      {
	writeError("No pointer size selected\n");
	return -3;
      }

    /* Parallelisation */
    if(parallel = lookupString(cgi,"parallel"))
      {
	if ( atoi(parallel) )
	  {
	    strcpy(executable+n,"_mp");
	  }
      }
    else
      {
	writeError("No parallelisation level set\n");
	return -4;
      }

    return 0;
  }

/* -- Job abort code --------------------------------------------------
*/

void printAbortOKButton( FILE *output, USER *user )
  {
    fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" JOB_CONTROL "\">\n" );
    fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->name );
    fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->password );
    fprintf( output, "<INPUT NAME=\"stage\" TYPE=\"hidden\" VALUE=\"11\">\n");
    fprintf( output, "<INPUT TYPE=\"submit\" VALUE=\"Abort Job\">\n" );
    fprintf( output, "</FORM>\n" );
  }


void printAbortAndClearButton( FILE *output, USER *user )
  {
    fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" JOB_CONTROL "\">\n" );
    fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->name );
    fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->password );
    fprintf( output, "<INPUT NAME=\"stage\" TYPE=\"hidden\" VALUE=\"12\">\n");
    fprintf( output, "<INPUT TYPE=\"submit\" VALUE=\"Abort Job/Clear Files\">\n" );
    fprintf( output, "</FORM>\n" );
  }


void printAbortCancelButton( FILE *output, USER *user )
  {
    fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" JOB_CONTROL "\">\n" );
    fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->name );
    fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->password );
    fprintf( output, "<INPUT TYPE=\"submit\" VALUE=\"Cancel\">\n" );
    fprintf( output, "</FORM>\n" );
  }


void writeAbortForm( CGI *cgi, FILE *output )
  {
    USER *user;
    JOB  *job;

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
	fprintf( stderr, "getJob() returned NULL in abortJob()\n" );
	disposeUser( user );
	return;
      }

    printJobHeading( output, user, job );
    fprintf( output, "<P>" );

    fprintf( output, "<FONT SIZE=\"+1\">\n" );
    if (strncmp(job->status,"running",4))
      {
	fprintf( output, "<CENTER>Are you sure you wish to abort the current job?</CENTER>\n" );
	fprintf( output, "</FONT><P>" );
	fprintf( output, "<CENTER>(The server doesn't think that it is running)</CENTER>\n" );
      }
    else
      {
	fprintf( output, "<CENTER>Are you sure you wish to abort the current job running on %s?"
		 "</CENTER>\n", job->host );
      }
    fprintf( output, "<P>\n" );

    fprintf( output, "<CENTER>" );
    fprintf( output, "<TABLE BORDER=\"0\" WIDTH=\"60%%\"><TR>" );
    fprintf( output, "<TD ALIGN=\"LEFT\">" );
    printAbortOKButton( output, user );
    fprintf( output, "</TD>" );
    fprintf( output, "<TD ALIGN=\"CENTER\">" );
    printAbortAndClearButton( output, user );
    fprintf( output, "</TD>" );
    fprintf( output, "<TD ALIGN=\"RIGHT\">" );
    printAbortCancelButton( output, user );
    fprintf( output, "</TD>" );
    fprintf( output, "</TR></TABLE>" );
    fprintf( output, "</CENTER>" );
  }

void abortJobCGI( CGI *cgi, FILE *output, HOSTLIST *hostList, int level )
  {
    USER *user;
    JOB  *job;
    HOST *host;
    HOSTITEM *hostItem;

    fprintf( output, "<H1>Abort Job</H1><P>" );
    fprintf( output, "<HR NOSHADE><P>" );

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
	writeError( "getJob() returned NULL in abortJobCGI()\n" );
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
	    fprintf( output, "<CENTER>\n" );
	    fprintf( output, "Flushing jobs/locks on %s.\n<p>", hostItem->host->hostname );
	    fprintf( output, "</CENTER>\n" );
	    
	    abortJob( output, hostItem->host, job, user, level);
	    hostItem = hostItem->nextHostItem;
	  }
      }
    else
      /* Abort on a named host */
      {
	/* rsh to hpc and do your business */
	if ( abortJob( output, host, job, user, level) )
	  {
	    fprintf( output, "<CENTER>\n" );
	    fprintf( output, "Error in abort job.\n" );
	    fprintf( output, "</CENTER>\n" );
	    fflush( output );
	  }
	else
	  {
	    fprintf( output, "<CENTER>\n" );
	    fprintf( output, "Your job has been aborted on %s.\n", host->hostname );
	    fprintf( output, "</CENTER>\n" );
	    fflush( output );
	  }
      }

    /* Bye! */
    fprintf( output, "<P>" );

    fprintf( output, "<CENTER>" );
    printOKButton( output, user );
    fprintf( output, "</CENTER>\n" );
  }


/* -- Job status code --------------------------------------------------
*/

void printStatusUpdateButton( FILE *output, USER *user )
  {
    fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" JOB_CONTROL "\">\n" );
    fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->name );
    fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->password );
    fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"stage\" VALUE=\"20\">\n");
    fprintf( output, "<INPUT TYPE=\"submit\" VALUE=\"Update\">\n" );
    fprintf( output, "</FORM>\n" );
  }


void printStatusReturnButton( FILE *output, USER *user )
  {
    fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" JOB_CONTROL "\">\n" );
    fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->name );
    fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->password );
    fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"stage\" VALUE=\"0\">\n");
    fprintf( output, "<INPUT TYPE=\"submit\" VALUE=\"Return\">\n" );
    fprintf( output, "</FORM>\n" );
  }


void writeStatusForm( CGI *cgi, FILE *output )
  {
    USER *user;
    JOB  *job;

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
	fprintf( stderr, "getJob() returned NULL in abortJob()\n" );
	disposeUser( user );
	return;
      }

    if (strncmp(job->status,"running",7) == 0)
      {
	borderGenerateHeader( output, user, True, JOB_CONTROL, "stage=20" );
      }
    else
      {
	borderGenerateHeader( output, user, False, JOB_CONTROL, "stage=20" );
      }

    printJobHeading( output, user, job );
    fprintf( output, "<P>" );

    fprintf( output, "<CENTER>" );
#ifndef NEW_CODE
    fprintf( output, "<TABLE BORDER=\"0\" WIDTH=\"40%%\">\n" );
    fprintf( output, "<TR><TD ALIGN=\"LEFT\" WIDTH=\"15%%\">\n" );
    printStatusUpdateButton( output, user );
    fprintf( output, "</TD>\n" );
    fprintf( output, "<TD ALIGN=\"RIGHT\" WIDTH=\"15%%\">\n" );
    printStatusReturnButton( output, user );
    fprintf( output, "</TR>\n</TABLE>\n" );
#else
    printStatusOKButton( output, user );
#endif
    fprintf( output, "</CENTER>" );
  }

/* -- initial menu --------------------------------------------------
*/

void writeJobControlForm( FILE *output, USER *user )
  {

    fprintf( output, "<H1>Job Control Options</H1><P>" );
    fprintf( output, "<HR NOSHADE><P>" );

    fprintf( output, "<FORM ACTION=\"" JOB_CONTROL "\" METHOD=\"POST\">\n" );
    fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">\n", user->name );
    fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"password\" VALUE=\"%s\">\n", user->password );
    fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"stage\" VALUE=\"1\">\n");
    fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"Submit job\">\n" );
    fprintf( output, "</FORM>\n" );
    fprintf( output, "<P>\n" );

    fprintf( output, "<FORM ACTION=\"" JOB_CONTROL "\" METHOD=\"POST\">\n" );
    fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">\n", user->name );
    fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"password\" VALUE=\"%s\">\n", user->password );
    fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"stage\" VALUE=\"10\">\n");
    fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"Abort job\">\n" );
    fprintf( output, "</FORM>\n" );
    fprintf( output, "<P>\n" );

    fprintf( output, "<FORM ACTION=\"" JOB_CONTROL "\" METHOD=\"POST\">\n" );
    fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">\n", user->name );
    fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"password\" VALUE=\"%s\">\n", user->password );
    fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"stage\" VALUE=\"20\">\n");
    fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"Job status\">\n" );
    fprintf( output, "</FORM>\n" );
    fprintf( output, "<P>\n" );

  }


/* -- Program Entry Point -----------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
    CGI *cgi;
    USER *user;
    char *passwordstring,*stagestring,*threadsstring,*userstring,
      *hostnamestring,*versionname;
    int stage;
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

    if(!(userstring = lookupString(cgi,"user")))
      {
	writeError("Invalid user");
	return 0;
      }

    if(!(passwordstring = lookupString(cgi,"password")))
      {
	writeError("Invalid password");
	return 0;
      }

    user = getUser( cgi );
    hostList = newHostList();
    readHostList( hostList, ESU_MASTER_HOST_LIST );

    if(!(stagestring = lookupString(cgi,"stage")))
      {
	stage=0;
      }
    else
      {
	stage=atoi(stagestring);
      }

    /* If anything has been canceled, then set to the default */
    if( lookupString( cgi, "cancel" ))
      {
	stage = 0;
      }


    /* Switch between the operation choices */
    switch (stage)
      {
      /* Main job control menu */
      case 0:
	{
	  borderGenerateHeader( stdout, user, False, NULL, NULL );
	  writeJobControlForm( stdout, user );
	  break;
	}

      /* Start of job submission -- select host */
      case 1:
	{ 
	  borderGenerateHeader( stdout, user, False, NULL, NULL );
	  writeHostForm( stdout, userstring, passwordstring, hostList );
	  break;
	}

      /* Stage 2 of job submission -- select process options */
      case 2:
	{
	  borderGenerateHeader( stdout, user, False, NULL, NULL );
	  hostnamestring = lookupString(cgi,"hostname");
	  if (hostnamestring == NULL)
	    {
	      writeError("Invalid hostname number");
	    }
	  else
	    {
	      writeVersionForm( stdout, userstring, passwordstring, hostnamestring, hostList );
	    }
	  break;
	}

      /* Stage 3 -- submit job */
      case 3:
	{
	  borderGenerateHeader( stdout, user, False, NULL, NULL );
	  if( getVersion( cgi, versionname, hostList ) == 0 )
	    {
	      hostnamestring = lookupString(cgi,"hostname");
	      threadsstring = lookupString(cgi,"threads");
	      if (hostnamestring == NULL)
		{
		  writeError("Invalid hostname number");
		}
	      else if (threadsstring == NULL)
		{
		  writeError("Invalid hostname number");
		}
	      else
		{
		  submitJobCGI( stdout, cgi , threadsstring, hostnamestring, versionname, hostList );
		}
	    }
	  break;
	}
	
      /* Stage 1 of job abortion -- form to stop job */
      case 10:
	{
	  borderGenerateHeader( stdout, user, True, NULL, NULL );
	  writeAbortForm( cgi, stdout );
	  break;
	}
	
      /* Stage 2 of job abortion -- kill the little fucker */
      case 11:
	{
	  borderGenerateHeader( stdout, user, False, NULL, NULL );
	  abortJobCGI( cgi, stdout, hostList, 0 );
	  break;
	}
	
      /* Stage 2a of job abortion -- kill proc, and blow away the files */
      case 12:
	{
	  borderGenerateHeader( stdout, user, False, NULL, NULL );
	  abortJobCGI( cgi, stdout, hostList, 1 );
	  break;
	}
	
	
      /* Job status */
      case 20:
	{
	  /* borderGenerateHeader( stdout, user, True, JOB_CONTROL, "stage=20" ); */
	  writeStatusForm( cgi, stdout );
	  break;
	}


      default:
	{
	  borderGenerateHeader( stdout, user, False, NULL, NULL );
	  writeError("Invalid stage number");
	  break;
	}
      } 
    
    borderGenerateFooter( stdout );
    fflush( stdout );

    return 0;
  }
