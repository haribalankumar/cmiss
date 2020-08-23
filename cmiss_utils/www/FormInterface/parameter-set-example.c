/*
 *  Parameter Set Example.c
 *
 *    Changes the users current set of files to a new example.
 */

/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stderr, etc */

#include <stdlib.h>
    /* For: malloc(), free(), NULL, etc */

#include <string.h>
    /* For: strlen(), strcpy(), strstr(), etc */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <grp.h>
    /* For: unlink() */

#include "directory.h"
    /* For: DIRECTORY, newDirectory(), disposeDirectory(), filterDirecory() */

#include "get-user.h"
    /* For: getUser() */

#include "user-managment.h"
    /* For: USER, disposeUser() */

#include "job.h"
    /* For: JOB getJob(), disposeJob() etc */

#include "cgi-decode.h"
    /* For: CGI handling functions */

#include "border.h"
    /* For: getBorderHeader() */

#include "transfer.h"
    /* For: copy() */

#include "webcmiss.h"
    /* For: file location defs */




char *generateSubName( char *dirname, char *subname )
  {
    char *fullName;
    int   length;

    /* Get the length of the fullName string */

    length = strlen( dirname ) + strlen( subname) ;

    /* alloc, and generate */
  
    fullName = malloc( 2 + length );

    strcpy( fullName, dirname );
    strcat( fullName, subname );

    return fullName;
  }



/* -- Module Methods ----------------------------------------------------------
*/

void copyFile( FILE *output, USER *user, char *sourceName, char *destinationName )
  {
    FILE    *source;
    FILE    *destination;
    char     buffer[512];
    char     command[512];
  
    source = fopen( sourceName, "rt" );
    if( source == NULL )
      {
	paramBorderWriteError( output, user, "Change Example: Error"
			       "Cannot access example directory" );
	exit(0);
      }

    destination = fopen( destinationName, "wt" );
    if( destination == NULL )
      {
	paramBorderWriteError( output, user, "Change Example: Error"
			       "Cannot access Input directory" );
	exit(0);
      }

    fgets( buffer, 512, source );
    while( !feof( source ) )
      {
	fputs( buffer, destination );

	fgets( buffer, 512, source );
      }

    fclose( source );
    fclose( destination );

    if(strstr(destinationName,".cgi"))
      {
	chmod( destinationName, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH );
	/*
	  sprintf(command,"chmod ugo+x %s",destinationName);
	  system(command);
	*/
      }
    if(strstr(destinationName,".gif"))
      {
	sprintf(command,"cp %s %s", sourceName, destinationName);
	system(command);
      }

    return;
  }


void copyFiles( FILE *output, USER *user, DIRECTORY *sourceDirectory, 
		char *destination, BOOLEAN test )
  {
    char         *fullSourceName;
    char         *fullDestinationName;
    FILEINFOITEM *fileInfoItem;
    FILEINFO     *fileInfo;
  
    fileInfoItem = sourceDirectory->fileInfoList->firstFileInfoItem;
    while( fileInfoItem != NULL )
      {
	fileInfo = fileInfoItem->fileInfo;

	if( fileInfo->type == NormalFileType )
	  {
	    /* skip known non-items */
	    if( strcmp( fileInfo->name, "description.html" ) != 0 &&      
		strcmp( fileInfo->name, "name.txt" ) != 0)
	      {
		fullSourceName = malloc( strlen( sourceDirectory->path ) + 
					 strlen( fileInfo->name ) + 2 );   /* + '/' + '\0' */
		if( fullSourceName == NULL )
		  {
		    paramBorderWriteError( output, user, "Change Example: Error" 
					   "Memory allocation failure in copyFiles() [1]" );
		    return;
		  }
		sprintf( fullSourceName, "%s/%s", sourceDirectory->path, 
			 fileInfo->name );
		
		fullDestinationName = malloc( strlen( destination ) + 
					      strlen( fileInfo->name ) + 2 );   /* + '/' + '\0' */
		if( fullDestinationName == NULL )
		  {
		    paramBorderWriteError( output, user, "Change Example: Error" 
					   "Memory allocation failure in copyFiles() [2]" );
		    free( fullSourceName );
		    return;
		  }
		sprintf( fullDestinationName, "%s/%s", destination, fileInfo->name );
		copy( fullSourceName, fullDestinationName );
		if (test) return;

		free( fullDestinationName );
		free( fullSourceName );
	      }
	  }

	fileInfoItem = fileInfoItem->nextFileInfoItem;
      }
    return;
  }


void removeFiles( FILE *output, USER *user, DIRECTORY *directory )
  {
    char         *fullFileName;
    FILEINFOITEM *fileInfoItem;
    FILEINFO     *fileInfo;
  
    fileInfoItem = directory->fileInfoList->firstFileInfoItem;
    while( fileInfoItem != NULL )
      {
	fileInfo = fileInfoItem->fileInfo;

	if( fileInfo->type == NormalFileType )
	  {
	    fullFileName = malloc( strlen( directory->path ) +
				   strlen( fileInfo->name ) + 2 );    /* + '/' + '\0' */
	    if( fullFileName == NULL )
	      {
		paramBorderWriteError( output, user, "Memory allocation failure in removeFiles()\n" );
		return;
	      }

	    sprintf( fullFileName, "%s/%s", directory->path, fileInfo->name );

	    remove( fullFileName );
	    free( fullFileName );
	  }

	fileInfoItem = fileInfoItem->nextFileInfoItem;
      }
  }


void paramSetExample( FILE *output, CGI *cgi, USER *user, JOB *job )
  {
    char      *example;
    char      *exampleDirectoryName;
    char      *buff;
    DIRECTORY *inputDirectory;
    DIRECTORY *exampleDirectory;
#ifdef SETGID
    struct group *cmiss;
#endif
    int len;

#ifdef SETGID
    /* Set GID to cmiss */
    if ((cmiss = getgrnam("cmiss")) == NULL)
      {
	paramBorderWriteError( output, user, "getgrnam() error in paramSetExample()\n" );
	return;
      }
    if (setgid(cmiss->gr_gid))
      {
	paramBorderWriteError( output, user, "setgid() error in paramSetExample()\n" );
	return;
      }
#endif

    /* Now set about copying examples */
    example = lookupString( cgi, "example" );
    if( example == NULL )
      {
	paramBorderWriteError( output, user, "No example name passed to paramSetExample()\n" );
	return;
      }

    /* Get the example directory */
    exampleDirectoryName = generateSubName( PARAM_DIR, example );

    if( exampleDirectoryName == NULL )
      {
	paramBorderWriteError( output, user, "generateExampleDirectoryNameNew() returned NULL in "
			       "paramSetExample()\n" );
	return;
    }

    exampleDirectory = newDirectory( exampleDirectoryName );
    if( exampleDirectory == NULL )
      {
	paramBorderWriteError( output, user, "newDirectory() returned NULL in paramSetExample()\n" );
	free( exampleDirectoryName );
	return;
      }

    /* Check we can copy the files */
    copyFiles( output, user, exampleDirectory, job->inputDirectoryName, True );

    /* Now open the destination (Input) directory, and clear it */
    inputDirectory = newDirectory( job->inputDirectoryName );
    if( inputDirectory == NULL )
      {
	paramBorderWriteError( output, user, "newDirectory() returned NULL in paramSetExample()\n" );
	free( exampleDirectoryName );
	return;
      }

    filterDirectory( inputDirectory, NULL );
    removeFiles( output, user, inputDirectory );
    disposeDirectory( inputDirectory );

    /* Copy the files across */
    filterDirectory( exampleDirectory, NULL );
    copyFiles( output, user, exampleDirectory, job->inputDirectoryName, False );

    disposeDirectory( exampleDirectory );

    /* Write out the job file */
    len = strlen(example) + 21;
    if ( ( job->description = malloc(len) ) == NULL )
      {
	paramBorderWriteError( output, user, "malloc fails in paramSetExample()\n" );
	return;
      }

    sprintf( job->description, "CMISS Example File %s", example );
    job->originalDirectoryName = exampleDirectoryName;

    writeJobFile( user, job );

    /* Forward onward to the example page */
    len = strlen(example) + 27;
    buff = malloc( len );
    if ( buff == NULL )
      {
	paramBorderWriteError( output, user, "malloc fails in paramSetExample()\n" );
	return;
      }

    sprintf( buff, "directory=input&example=%s", example );

    forwardUserTag( PARAM_EXAMPLE, user, buff );

    free( buff );
  }




/* -- Program Entry Point -----------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
    USER *user;
    JOB   *job;
    CGI   *cgi;

    cgi = getCGIEnvironment( argc, argv );
    if( cgi == NULL )
      {
	printf( "<HR><P><B>ERROR:</B> getCGIEnvironment() failed in main()</P>"
		"<P>Email "WEBMASTER"<I>immediately</I></P>\n" );
	return -1;
      }

    user = getUser( cgi );
    if( user == NULL )
      {
	/* Not always an error, the user could simply have supplied an incorrect
	   password. getUser() will report back those conditions to the user
	   (one hopes...) */
	return -1;
      }

    job = getJob( user );
    if( job == NULL )
      {
	fprintf( stderr, "getJob() returned NULL in paramSetExample()\n" );
	disposeUser( user );
	return -1;
      }

    paramSetExample( stdout, cgi, user, job );

    return 0;
  }

