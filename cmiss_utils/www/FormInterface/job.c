/*
 *  Job.c
 *
 *    The JOB datatype.
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stdout, FILE, etc */

#include <stdlib.h>
    /* For: malloc(), free(), NULL */

#include <sys/types.h>
#include <sys/stat.h>
    /* For: stat() */

#include <string.h>
    /* For: strlen(), strchr(), strstr(), etc */

#include "get-user.h"
    /* For: getUser() */

#include "cgi-decode.h"
    /* For: CGI handling functions */

#include "directory.h"
    /* For: DIRECTORY, FILEINFO, etc */

#include "groups.h"
    /* For: GROUP, newGroup(), getTag(), etc */
#include "job.h"
    /* For: Our public interface */

#include "webcmiss.h"

/* -- Private Module Methods --------------------------------------------------
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


static void stripTrailingSpace( char *string )
  {
  char *p;
  char *lastSpace;

  p = string + strlen( string );
  lastSpace = p;

  while( p != string )
    {
    if( *p == '\r' || *p == '\n' || *p == '\t' || *p == ' ' || *p == '\0' )
      lastSpace = p;
    else
      break;

    p--;
    }

  *lastSpace = '\0';
  }



/* -- Directory Filter Methods ------------------------------------------------
*/

void filterDirectory( DIRECTORY *directory, void *filterPatterns )
  {
  FILEINFOITEM  *fileInfoItem;
  FILEINFOITEM **fileInfoItemPtr;
  FILEINFO      *fileInfo;

	if (filterPatterns != NULL) {
	}

  fileInfoItemPtr = &(directory->fileInfoList->firstFileInfoItem);
  fileInfoItem    = *fileInfoItemPtr;
  while( fileInfoItem != NULL )
    {
    fileInfo = fileInfoItem->fileInfo;

    /* simplistic filter, make a RE comparison later */
    if( fileInfo->name[0] == '.' )
      {
      *fileInfoItemPtr = fileInfoItem->nextFileInfoItem;
      disposeFileInfoItem( fileInfoItem );
      fileInfoItem = *fileInfoItemPtr;
      }
    else
      {
      fileInfoItemPtr = &(fileInfoItem->nextFileInfoItem);
      fileInfoItem    = *fileInfoItemPtr;
      }
    }
  }


void sortDirectory( DIRECTORY *directory )
  {
		if (directory != NULL) {
		}
  }


/* -- Public Methods ----------------------------------------------------------
*/

JOB *newJob( void )
  {
  JOB *tempJob;

  tempJob = malloc( sizeof( JOB ) );
  if( tempJob == NULL )
    {
    fprintf( stderr, "memory allocation failure in newJob()\n" );
    return NULL;
    }

  tempJob->description = NULL;
  tempJob->originalDirectoryName = NULL;
  tempJob->status = NULL;
  tempJob->host = NULL;
  tempJob->pid = NULL;
  tempJob->startTime = NULL;
  tempJob->inputDirectoryName = NULL;
  tempJob->outputDirectoryName = NULL;
  tempJob->inputFiles = NULL;
  tempJob->outputFiles = NULL;

  return tempJob;
  }


void disposeJob( JOB *job )
  {
  if( job->description != NULL )
    free( job->description );

  if( job->originalDirectoryName != NULL )
    free( job->originalDirectoryName );

  if( job->status != NULL )
    free( job->status );

  if( job->host != NULL )
    free( job->host );

  if( job->pid != NULL )
    free( job->pid );

  if( job->startTime != NULL )
    free( job->startTime );

  if( job->inputDirectoryName != NULL )
    free( job->inputDirectoryName );

  if( job->outputDirectoryName != NULL )
    free( job->outputDirectoryName );

  if( job->inputFiles != NULL )
    disposeDirectory( job->inputFiles );

  if( job->outputFiles != NULL )
    disposeDirectory( job->outputFiles );

  free( job );
  }


JOB *getJob( USER *user )
  { 
  char      *fileName;
  GROUPFILE *file;
  RULES     *rules;
  GROUP     *group;
  TAG       *tag;
  char      *jobDescription;
  char      *originalDirectoryName;
  char      *status;
  char      *host;
  char      *pid;
  char      *startTime;
  char      *inputDirectoryName;
  char      *outputDirectoryName;

  char      *directoryName;
  DIRECTORY *directory;
  JOB       *job;

  /* Constants */
  /*Now defined in Makefile*/
/*  const char *RuleFileName = "/usr/people/poor/CMISS/FormInterface/job.rules";*/
/*  const char *JobFileName  = "job.info";*/

  job = newJob();
  if( job == NULL )
    {
    fprintf( stderr, "newJob() returned NULL in getJob()\n" );
    return NULL;
    }


  /* -- Read the job.info file --------------------------------------------- */

  fileName = malloc( strlen( user->homeDirectory ) + strlen( JOB_FILE_NAME ) +
    2 );   /*  "/user/bob" + "/" + "job.info" + "\0" */
  if( fileName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in getJob()\n" );
    disposeJob( job );
    return NULL;
    }

  sprintf( fileName, "%s/%s", user->homeDirectory, JOB_FILE_NAME );

  file = newGroupFile( fileName );
  if( file == NULL )
    {
    fprintf( stderr, "newGroupFile() returned NULL in getJob()\n" );
    free( fileName );
    disposeJob( job );
    return NULL;
    }

  rules = newRules();
  if( rules == NULL )
    {
    fprintf( stderr, "newRules() returned NULL in getJob()\n" );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }
  
  parseRuleFile( rules, ESU_RULE_FILE_NAME );
    /* -- NB. No error checking is preformed at this point */

  group = getGroup( file, rules );
  if( group == NULL )
    {
    fprintf( stderr, "getGroup() returned NULL in getJob()\n" );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }

  tag = lookupTag( group, "job" );
  if( tag == NULL )
    {
    fprintf( stderr, "getTag() returned NULL in getJob()\n" );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }

  jobDescription = newStringCopy( tag->body );
  if( jobDescription == NULL )
    {
    fprintf( stderr, "newStringCopy() returned NULL in getJob()\n" );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }

  tag = lookupTag( group, "original-files" );
  if( tag == NULL )
    {
    fprintf( stderr, "getTag() returned NULL in getJob()\n" );
    free( jobDescription );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }

  originalDirectoryName = newStringCopy( tag->body );
  if( originalDirectoryName == NULL )
    {
    fprintf( stderr, "newStringCopy() returned NULL in getJob()\n" );
    free( jobDescription );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }

  tag = lookupTag( group, "status" );
  if( tag == NULL )
    {
    fprintf( stderr, "getTag() returned NULL in getJob()\n" );
    free( originalDirectoryName );
    free( jobDescription );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }

  status = newStringCopy( tag->body );
  if( status == NULL )
    {
    fprintf( stderr, "newStringCopy() returned NULL in getJob()\n" );
    free( originalDirectoryName );
    free( jobDescription );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }

  tag = lookupTag( group, "host" );
  if( tag == NULL )
    {
    fprintf( stderr, "getTag() returned NULL in getJob()\n" );
    free( status );
    free( originalDirectoryName );
    free( jobDescription );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }

  host = newStringCopy( tag->body );
  if( status == NULL )
    {
    fprintf( stderr, "newStringCopy() returned NULL in getJob()\n" );
    free( status );
    free( originalDirectoryName );
    free( jobDescription );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }

  tag = lookupTag( group, "pid" );
  if( tag == NULL )
    {
    fprintf( stderr, "getTag() returned NULL in getJob()\n" );
    free( host );
    free( status );
    free( originalDirectoryName );
    free( jobDescription );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }

  pid = newStringCopy( tag->body );
  if( pid == NULL )
    {
    fprintf( stderr, "newStringCopy() returned NULL in getJob()\n" );
    free( host );
    free( status );
    free( originalDirectoryName );
    free( jobDescription );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }

  tag = lookupTag( group, "startTime" );
  if( tag == NULL )
    {
    fprintf( stderr, "getTag() returned NULL in getJob()\n" );
    free( pid );
    free( host );
    free( status );
    free( originalDirectoryName );
    free( jobDescription );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }

  startTime = newStringCopy( tag->body );
  if( startTime == NULL )
    {
    fprintf( stderr, "newStringCopy() returned NULL in getJob()\n" );
    free( pid );
    free( host );
    free( status );
    free( originalDirectoryName );
    free( jobDescription );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }



  stripTrailingSpace( originalDirectoryName );

  inputDirectoryName = malloc( strlen( user->homeDirectory ) + 7 );
    /* + "/Input" + '\0' */
  if( inputDirectoryName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in getJob()\n" );
    fprintf( stderr, "newStringCopy() returned NULL in getJob()\n" );
    free( startTime );
    free( pid );
    free( host );
    free( status );
    free( originalDirectoryName );
    free( jobDescription );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }
  sprintf( inputDirectoryName, "%s/Input", user->homeDirectory );

  outputDirectoryName = malloc( strlen( user->homeDirectory ) + 8 );
    /* + "/Output" + '\0' */
  if( inputDirectoryName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in getJob()\n" );
    fprintf( stderr, "newStringCopy() returned NULL in getJob()\n" );
    free( startTime );
    free( pid );
    free( host );
    free( status );
    free( originalDirectoryName );
    free( inputDirectoryName );
    free( jobDescription );
    free( fileName );
    disposeGroupFile( file );
    disposeJob( job );
    return NULL;
    }
  sprintf( outputDirectoryName, "%s/Output", user->homeDirectory );

  job->description = jobDescription;
  job->originalDirectoryName = originalDirectoryName;
  job->inputDirectoryName = inputDirectoryName;
  job->outputDirectoryName = outputDirectoryName;
  job->status = status;
  job->host = host;
  job->pid = pid;
  job->startTime = startTime;

  free( fileName );
  disposeRules( rules );
  disposeGroupFile( file );


  /* -- Read the Directories ----------------------------------------------- */

  directoryName = malloc( strlen( user->homeDirectory ) + 8 );
                                           /* + '/Output' + '\0' */
  if( directoryName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in getJob()\n" );
    disposeJob( job );
    return NULL;
    }

  /* Input Files */

  sprintf( directoryName, "%s/Input", user->homeDirectory );
  directory = newDirectory( directoryName );
  if( directory == NULL )  
    {
    printf( "newDirectory() returned NULL in getJob() [1]\n" );
    free( directoryName );
    disposeJob( job );
    return NULL;
    }

  filterDirectory( directory, NULL );
  job->inputFiles = directory;
  
  /* Output Files */
  sprintf( directoryName, "%s/Output", user->homeDirectory );
  directory = newDirectory( directoryName );
  if( directory == NULL )  
    {
    printf( "newDirectory() returned NULL in getJob() [2]\n" );
    free( directoryName );
    disposeDirectory( job->inputFiles );
    disposeJob( job );
    return NULL;
    }

  filterDirectory( directory, NULL );
  job->outputFiles = directory;

  return job;
  }


void writeJobFile( USER *user, JOB *job )
  {
  FILE *jobFile;
  char *jobFileName;
	struct stat fileStat;

  jobFileName = malloc( strlen( user->homeDirectory ) + 10 );
      /* + "job.info" + '/' + '\0' */
  if( jobFileName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in writeJobFile()\n" );
    return;
    }

  sprintf( jobFileName, "%s/job.info", user->homeDirectory );
	/* If the job file exists, remove it */
  if (stat( jobFileName, &fileStat) == 0) {
		remove( jobFileName );
	}

  jobFile = fopen( jobFileName, "wt" );
  if( jobFile == NULL )
    {
    fprintf( stderr, "fopen() returned NULL in writeJobFile()\n" );
    free( jobFileName );
    return;
    }

	stripTrailingSpace(job->description);
	stripTrailingSpace(job->originalDirectoryName);
	stripTrailingSpace(job->status);

  fprintf( jobFile, "Job: %s\n", job->description );
  fprintf( jobFile, "Original-Files: %s\n", job->originalDirectoryName );
	fprintf( jobFile, "Status: %s\n", job->status );
	if (strncmp(job->status,"running",7) == 0)
		{
		stripTrailingSpace(job->host);
		stripTrailingSpace(job->pid);
		stripTrailingSpace(job->startTime);

    fprintf( jobFile, "Host: %s\n", job->host );
    fprintf( jobFile, "Pid: %s\n", job->pid );
    fprintf( jobFile, "StartTime: %s\n", job->startTime );
		}
	else
		{
    fprintf( jobFile, "Host: -\n" );
    fprintf( jobFile, "Pid: 0\n" );
    fprintf( jobFile, "StartTime: 0\n" );
		}

  fclose( jobFile );
  }


void printJobHeading( FILE *output, USER *user, JOB *job )
  {
	time_t currentTime,startingTime;
	char *time_p;

  fprintf( output, "<CENTER>" );
  fprintf( output, "<TABLE BGCOLOR=\"#ddbb88\" BORDER=\"0\" "
    "CELLPADDING=\"3\" CELLSPACING=\"0\"><TR><TD>" );
  fprintf( output, "<TABLE BGCOLOR=\"#ffffcc\" BORDER=\"0\" "
    "CELLPADDING=\"3\" CELLSPACING=\"0\">" );
  fprintf( output, "<TR><TD ALIGN=\"RIGHT\"><FONT SIZE=\"+1\" FACE=\"Helvetica\">"
					"<B>Username:</B></FONT></TD><TD><FONT SIZE=\"+1\" FACE=\"Helvetica\">%s</FONT></TD></TR>",
					user->name );
  fprintf( output, "<TR><TD ALIGN=\"RIGHT\"><FONT SIZE=\"+1\" FACE=\"Helvetica\">"
					"<B>Current Job:</FONT></B></TD><TD><FONT SIZE=\"+1\" FACE=\"Helvetica\">%s</FONT></TD></TR>",
					job->description );
	if (strncmp(job->status,"running",7) == 0)
		{
	  currentTime = time(NULL);
		startingTime = atoi(job->startTime);
		time_p = ctime(&startingTime);
		time_p[24] = '\0';

		fprintf( output, "<TR><TD ALIGN=\"RIGHT\"><FONT SIZE=\"+1\" FACE=\"Helvetica\">"
						"<B>Status:</FONT></B></TD><TD><FONT SIZE=\"+1\" FACE=\"Helvetica\">%s</FONT></TD></TR>",
						"Running" );
		fprintf( output, "<TR><TD ALIGN=\"RIGHT\"><FONT SIZE=\"+1\" FACE=\"Helvetica\">"
						"<B>Host:</FONT></B></TD><TD><FONT SIZE=\"+1\" FACE=\"Helvetica\">%s</FONT></TD></TR>",
						job->host );
		fprintf( output, "<TR><TD ALIGN=\"RIGHT\"><FONT SIZE=\"+1\" FACE=\"Helvetica\">"
						"<B>PID:</FONT></B></TD><TD><FONT SIZE=\"+1\" FACE=\"Helvetica\">%s</FONT></TD></TR>",
						job->pid );
		fprintf( output, "<TR><TD ALIGN=\"RIGHT\"><FONT SIZE=\"+1\" FACE=\"Helvetica\">"
						"<B>Start Time:</FONT></B></TD><TD><FONT SIZE=\"+1\" FACE=\"Helvetica\">%s</FONT></TD></TR>",
						time_p );
		fprintf( output, "<TR><TD ALIGN=\"RIGHT\"><FONT SIZE=\"+1\" FACE=\"Helvetica\">"
						"<B>Running Time:</FONT></B></TD><TD><FONT SIZE=\"+1\" FACE=\"Helvetica\">%d seconds</FONT></TD></TR>",
						currentTime - startingTime );
	  }
	else
		{
		fprintf( output, "<TR><TD ALIGN=\"RIGHT\"><FONT SIZE=\"+1\" FACE=\"Helvetica\">"
						"<B>Status:</FONT></B></TD><TD><FONT SIZE=\"+1\" FACE=\"Helvetica\">%s</FONT></TD></TR>",
						"Not-Running" );
	  }
  fprintf( output, "</TABLE>" );
  fprintf( output, "</TD></TR></TABLE>" );
  fprintf( output, "</CENTER>" );

	fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">",user->name);
	fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">",user->password);
  }



void printBriefJobHeading( FILE *output, JOB *job )
{
  time_t currentTime,startingTime;
  char *time_p;

  fprintf( output, "<CENTER>" );
  fprintf( output, "<TABLE BGCOLOR=\"#ddbb88\" BORDER=\"0\" "
           "CELLPADDING=\"3\" CELLSPACING=\"0\"><TR><TD>" );
  fprintf( output, "<TABLE BGCOLOR=\"#ffffcc\" BORDER=\"0\" "
           "CELLPADDING=\"3\" CELLSPACING=\"0\">" );

  if (strncmp(job->status,"running",7) == 0) {
    currentTime = time(NULL);
    startingTime = atoi(job->startTime);
    time_p = ctime(&startingTime);
    time_p[24] = '\0';

    fprintf( output, "<TR><TD ALIGN=\"RIGHT\"><FONT SIZE=\"+1\" FACE=\"Helvetica\">"
             "<B>Status:</FONT></B></TD><TD><FONT SIZE=\"+1\" FACE=\"Helvetica\">%s</FONT></TD></TR>",
             "Running" );
    fprintf( output, "<TR><TD ALIGN=\"RIGHT\"><FONT SIZE=\"+1\" FACE=\"Helvetica\">"
             "<B>Running Time:</FONT></B></TD><TD><FONT SIZE=\"+1\" FACE=\"Helvetica\">%d seconds</FONT></TD></TR>",
             currentTime - startingTime );
  }
  else {
    fprintf( output, "<TR><TD ALIGN=\"RIGHT\"><FONT SIZE=\"+1\" FACE=\"Helvetica\">"
				     "<B>Status:</FONT></B></TD><TD><FONT SIZE=\"+1\" FACE=\"Helvetica\">%s</FONT></TD></TR>",
					   "Not-Running" );
  }
  
  fprintf( output, "</TABLE>" );
  fprintf( output, "</TD></TR></TABLE>" );
  fprintf( output, "</CENTER>" );
}




