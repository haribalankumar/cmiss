/*
 *  Change Example.c
 *
 *    Changes the users current st of files to a new set from the
 *    examples directories.
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


/*
 *  The following code assumes that the name passed in will be in the
 *  form "2a3", which means that the resulting example file will
 *  be found in ".../example_2/example_2a/example_2a3/example_2a3.com"
 *
 *  This is a bit messy, but never mind.
 */
char *generateExampleDirectoryName( char *examples_path, char *exampleName )
  {
  char *fullName;
  char *index;
  int   length;
  int   fullLength;
  int   i;

  /* Handle backwards compatibility (welcome to the rest of your life) */
  index = strstr( exampleName, "example_" );
  if( index != NULL )
    {
    exampleName = index + strlen( "example_" ); 

    while( *index != '\0' && *index != '.' ) 
      index++;

    *index = '\0'; 
    }

  length = strlen( exampleName );
  fullLength = 9 * ( length + 1 ) + length * length; /* Over-allocates */
    /* Actually the previous line over-allocated before, but now that
       I've commented out the strcat line below, it vastly over-allocates.
       This is in need of a fix */

  fullName = malloc( strlen( examples_path ) + 1 + fullLength );

  strcpy( fullName, examples_path );

  for( i = 0; i < length; i++ )
    {
    /*strcat( fullName, "example_" );*/
    strncat( fullName, exampleName, i+1 );
    strcat( fullName, "/" );
    }

  /*
  strcat( fullName, "example_" );
  strcat( fullName, exampleName );
  strcat( fullName, ".com" );
  */

  return fullName;
  }


/*
 *  New version of generateExampleDirectoryName. The name passed in will
 *  be in the form "2a3", and the resulting example files will be found
 *  in the ".../2/2a/2a3/cmiss_input/" directory.
 */
char *generateExampleDirectoryNameNew( char *examples_path, char *exampleName )
  {
  char *fullName;
  int   length;
  int   fullLength;
  int   i;
  int   len;


  /* Get the length of the fullName string */

  length = strlen( exampleName );

  len = 0;
  for ( i = 0; i < length ; i++ )
    {
    len += (i+2);
    }
  
  fullLength = len + 14;


  /* alloc, and generate */
  
  fullName = malloc( strlen( examples_path ) + 1 + fullLength );

  strcpy( fullName, examples_path );
  strcat( fullName, "/" );

  for( i = 0; i < length; i++ )
    {
    strncat( fullName, exampleName, i+1 );
    strcat( fullName, "/" );
    }

  strcat( fullName, "cmiss_input/" );

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
		borderWriteError( output, user, "Change Example: Error"
										 "Cannot access example directory" );
    exit(0);
    }

  destination = fopen( destinationName, "wt" );
  if( destination == NULL )
    {
		borderWriteError( output, user, "Change Example: Error"
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
					borderWriteError( output, user, "Change Example: Error" 
													 "Memory allocation failure in copyFiles() [1]" );
          return;
          }
        sprintf( fullSourceName, "%s/%s", sourceDirectory->path, 
          fileInfo->name );

        fullDestinationName = malloc( strlen( destination ) + 
          strlen( fileInfo->name ) + 2 );   /* + '/' + '\0' */
        if( fullDestinationName == NULL )
          {
					borderWriteError( output, user, "Change Example: Error" 
													 "Memory allocation failure in copyFiles() [2]" );
          free( fullSourceName );
          return;
          }
        sprintf( fullDestinationName, "%s/%s", destination, fileInfo->name );
#ifdef OLD_CODE
        copyFile( output, user, fullSourceName, fullDestinationName );
#else
        copy( fullSourceName, fullDestinationName );
#endif
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
				borderWriteError( output, user, "Memory allocation failure in removeFiles()\n" );
        return;
        }

      sprintf( fullFileName, "%s/%s", directory->path, fileInfo->name );

      remove( fullFileName );
      free( fullFileName );
      }

    fileInfoItem = fileInfoItem->nextFileInfoItem;
    }
  }


void changeExample( FILE *output, CGI *cgi, USER *user, JOB *job )
  {
  char      Input_dir[512];
  char      *exampleNumber;
  char      *exampleDirectoryName;
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
		borderWriteError( output, user, "getgrnam() error in changeExample()\n" );
    return;
		}
	if (setgid(cmiss->gr_gid))
		{
		borderWriteError( output, user, "setgid() error in changeExample()\n" );
    return;
		}
#endif

	/* Now set about copying examples */
  exampleNumber = lookupString( cgi, "example" );
  if( exampleNumber == NULL )
    {
		borderWriteError( output, user, "No example number passed to changeExample()\n" );
    return;
    }

	/* Get the example directory */
  exampleDirectoryName = generateExampleDirectoryNameNew( ESU_EXAMPLES_PATH,
    exampleNumber );
  if( exampleDirectoryName == NULL )
    {
		borderWriteError( output, user, "generateExampleDirectoryNameNew() returned NULL in "
      "changeExample()\n" );
    return;
    }

  exampleDirectory = newDirectory( exampleDirectoryName );
  if( exampleDirectory == NULL )
    {
		borderWriteError( output, user, "newDirectory() returned NULL in changeExample()\n" );
    free( exampleDirectoryName );
    return;
    }
	/* Check we can copy the files */
  copyFiles( output, user, exampleDirectory, job->inputDirectoryName, True );

	/* Now open the destination (Input) directory, and clear it */
  inputDirectory = newDirectory( job->inputDirectoryName );
  if( inputDirectory == NULL )
    {
		borderWriteError( output, user, "newDirectory() returned NULL in changeExample()\n" );
    free( exampleDirectoryName );
    return;
    }

  filterDirectory( inputDirectory, NULL );
  removeFiles( output, user, inputDirectory );
  disposeDirectory( inputDirectory );

	/* Copy the files across */
  filterDirectory( exampleDirectory, NULL );
  copyFiles( output, user, exampleDirectory, job->inputDirectoryName, False );

  strncpy(Input_dir,job->inputDirectoryName,strlen(job->inputDirectoryName)-5);

  disposeDirectory( exampleDirectory );

	/* Write out the job file */
	len = strlen(exampleNumber) + 21;
	if ( ( job->description = malloc(len) ) == NULL )
		{
		borderWriteError( output, user, "malloc fails in changeExample()\n" );
    return;
		}

	sprintf( job->description, "CMISS Example File %s", exampleNumber );
	job->originalDirectoryName = exampleDirectoryName;

	writeJobFile( user, job );

	/* Forward onward to the file lister */
	/* forwardUser( LIST_PROBLEM, user ); */
	forwardUserDir( LIST_PROBLEM, user, "input" );
  }


/* -- Change Form ----------------------------------------------------------
*/

void printChangeForm( FILE *output, USER *user, JOB *job )
  {
  fprintf( stdout, "<HR NOSHADE><P>\n" );
  printJobHeading( output, user, job );
  fprintf( stdout, "<P>\n" );

  fprintf( output, "<FONT SIZE=\"+1\">\n" );
	fprintf( output, "Please choose a new example number (e.g. 144, 148, 336, d21, h2, h3, b152, b221). \n" );
	fprintf( output, "See the <A HREF=\"" EXAMPLE_FILES_URL "\"> CMISS online example files</A> \n" );
	fprintf( output, "for a list of files to choose from.\n" );
  fprintf( output, "</FONT><P>\n" );

  fprintf( output, "<CENTER>\n" );
  fprintf( output, "<TABLE BORDER=\"0\" WIDTH=\"80%%\"><TR>\n" );
  fprintf( output, "<TD ALIGN=\"CENTER\" COLSPAN=\"2\">\n" );

  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" CHANGE_EXAMPLE "\">\n" );
  fprintf( output, "<FONT SIZE=\"+1\">Please input new example number:</FONT> \n" );
  fprintf( output, "<INPUT NAME=\"example\" TYPE=\"text\" VALUE=\"\" SIZE=\"12\">\n" );
  fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->name );
  fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->password );
  fprintf( output, "<INPUT NAME=\"stage\" TYPE=\"hidden\" VALUE=\"1\">\n" );
  fprintf( output, "</TD></TR>\n" );

  fprintf( output, "<TR>\n" );
  fprintf( output, "<TD ALIGN=\"LEFT\">\n" );
  fprintf( output, "<INPUT TYPE=\"submit\" VALUE=\"Change Example\">\n" );
  fprintf( output, "</FORM>\n" );
  fprintf( output, "</TD>\n" );

  fprintf( output, "<TD ALIGN=\"RIGHT\">\n" );
  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" LIST_PROBLEM "\">\n" );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">\n", user->name );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"password\" VALUE=\"%s\">\n", user->password );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"directory\" VALUE=\"input\">\n" );
  fprintf( output, "<INPUT TYPE=\"submit\" VALUE=\"Cancel\">\n" );
  fprintf( output, "</FORM>\n" );

  fprintf( output, "</TD>\n" );  
  fprintf( output, "</TR>\n" );
  fprintf( output, "</TABLE>\n" );
  fprintf( output, "</CENTER>\n" );

  fprintf( output, "<FONT SIZE=\"+1\" COLOR=\"#990000\">Note Well:</FONT>\n" );
	fprintf( output, "<FONT SIZE=\"+1\">\n" );
  fprintf( output, "Changing to another example will <B>permanantly overwrite</B> your \n" );
  fprintf( output, "existing files and will lose any changes you have made to them.\n" );
  fprintf( output, "</FONT><P>\n" );
  }


/* -- Program Entry Point -----------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
	USER *user;
	JOB   *job;
  CGI   *cgi;
	char *stagestring;
	int  stage;

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
    fprintf( stderr, "getJob() returned NULL in changeExample()\n" );
    disposeUser( user );
    return -1;
    }

  if(!(stagestring = lookupString(cgi,"stage")))
    {
    stage=0;
    }
  else
    {
    stage=atoi(stagestring);
    }


	/* Switch between the operation stages */
  switch (stage)
    {
		/* Example change menu */
  	case 0:
			{
			borderGenerateHeader( stdout, user, False, NULL, NULL );
			printChangeForm( stdout, user, job );
			borderGenerateFooter( stdout );
			break;
	    }

		/* Switch examples */
    case 1:
      {
			changeExample( stdout, cgi, user, job );
      break;
      }


    default:
      {
			borderGenerateHeader( stdout, user, False, NULL, NULL );
      writeError("Invalid stage number");
			borderGenerateFooter( stdout );
      break;
      }
    } 
  

  return 0;
  }

