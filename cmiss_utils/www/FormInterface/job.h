/*
 *  Job.h
 *
 *    The JOB datatype's public interface.
 */

#ifndef JOB_H
#define JOB_H


/* -- Required Include Files --------------------------------------------------
*/

#include <stdio.h>
    /* For: FILE */

#include "user-managment.h" 
    /* For: USER */

#include "directory.h"
    /* For: DIRECTORY */

#include "groups.h"
    /* For: GROUP, newGroup(), getTag(), etc */


/* -- Module Datatypes  -------------------------------------------------------
*/

typedef struct JOB
  {
  char      *description;
  char      *originalDirectoryName;
	char      *status;
	char      *host;
	char      *pid;
	char      *startTime;
  char      *inputDirectoryName;
  char      *outputDirectoryName;
  DIRECTORY *inputFiles;
  DIRECTORY *outputFiles;
  }
JOB;



/* -- Public Module Methods ---------------------------------------------------
*/

JOB *newJob( void );
void disposeJob( JOB *job );

JOB *getJob( USER *user );

void writeJobFile( USER *user, JOB *job );
void printJobHeading( FILE *output, USER *user, JOB *job );
void printBriefJobHeading( FILE *output, JOB *job );


#endif
