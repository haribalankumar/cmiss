/*
 *  Submission test.c
 *  
 *    Code to debug the submission server.
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), fflush(), stdout, stderr, etc */

#include <stdlib.h>
    /* For: malloc(), free(), NULL */

#include <string.h>
    /* For: strlen(), strcpy() */

#include <unistd.h>
    /* For: chdir() */

#include "user-managment.h"
    /* For: USER */

#include "get-user.h"
    /* For: getUser() */

#include "job.h"
    /* For: JOB */

#include "host.h"
    /* For: HOST, HOSTITEM etc */

#include "transfer.h"
    /* For: Socket communication */


/* -- Module Constants --------------------------------------------------------
*/
#include "webcmiss.h"
#define BUFFER_SIZE     1024


void abortJobTest( FILE *output, HOST *host, USER *user, JOB *job )
{
	/* rsh to hpc and do your business */
	if ( abortJob( output, host, job, user) ) {
	  fprintf( output, "abortJob() aborted in abortJobCGI()\n" );
	}
}

void submitJobTest( FILE *output, HOST *host, USER *user, JOB *job, char *threadstring, 
									  char *versionname, char *indir  )
{
	int ierr;

  if (chdir(indir)) {
    fprintf( output, "submitJobCGI: cannot cd to \"%s\"\n",indir );
    return;
  }

	/* Submit the job */
	ierr = submitJob( output, versionname, threadstring, host, job, user);

	switch(ierr)
  {
    case 0:
    {
      fprintf( output, "Your job has been submitted to %s.\n", host->hostname );
      break;
    }
    case 1:
    {
      fprintf( output, "A job is already running on %s.\n", host->hostname );
      return;
    }
    default:
    {
      fprintf( output, "Unknown return code from submitJob(): %d\n", ierr );
			fflush( output );
      break;
    }
  }
}


main(int argc, char **argv)
{
	HOSTLIST *hostList;
	HOST *host;
	USERLIST *userList;
	USER *user;

	JOB  *job;

	char *threadstring = "1";
	char *versionname  = "cm_32";
	char *hostBaseName = ESU_MASTER_HOST_LIST;
	char *passBaseName = ESU_MASTER_USER_LIST;
	char *indir        = ESU_WEB "people/norris/Input/";

	char *username = "norris";
	char *hostname = "hpc2";

	int i, abort = 0;

	for (i=1; i<argc; i++) {
		if (strncmp(argv[i],"-host",5) == 0) {
			i++;
      if (argv[i] != NULL) {
				hostname = argv[i];
			}
		}
		else if (strncmp(argv[i],"-user",5) == 0) {
			i++;
      if (argv[i] != NULL) {
				username = argv[i];
			}
		}
		else if (strncmp(argv[i],"-exe",4) == 0) {
			i++;
      if (argv[i] != NULL) {
				versionname = argv[i];
			}
		}

		else if (strncmp(argv[i],"-abort",6) == 0) {
			abort = 1;
		}
		else if (strncmp(argv[i],"-submit",7) == 0) {
			abort = 0;
		}
		else {
			fprintf(stderr, "Unknown flag \"%s\"\n", argv[i]);
			exit(-1);
		}
	}


	/* Get info from databases */
	hostList = newHostList();
	readHostList( hostList, hostBaseName );
	host = lookupHost( hostname, hostList );

	userList = newUserList();
	readUserList( userList, passBaseName );
	user = lookupUser( username, userList );

	job = newJob();
  job->description              = "test run";
  job->originalDirectoryName    = "";
	job->status                   = NULL;
	job->host                     = NULL;
	job->pid                      = NULL;
	job->startTime                = NULL;
  job->inputDirectoryName       = ESU_WEB "people/norris/Input";
  job->outputDirectoryName      = ESU_WEB "people/norris/Output";


	if (abort) {
		abortJobTest( stdout, host, user, job );
	}
	else {
		submitJobTest( stdout, host, user, job, threadstring, versionname, indir );
	}

	return 0;
}
