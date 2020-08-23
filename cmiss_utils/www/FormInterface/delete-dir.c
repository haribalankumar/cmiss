/*
 *  Delete Dir.c
 *
 *    Asks whether to delete all the files in a directory.
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), etderr, etc */

#include <stdlib.h>
    /* For: malloc() free(), NULL */

#include <string.h>
    /* For: strlen(), strchr(), strstr(), etc */

#include "get-user.h"
    /* For: getUser() */

#include "job.h"
    /* For: JOB, getJob(), etc */

#include "cgi-decode.h"
    /* For: CGI handling functions */

#include "directory.h"
    /* For: DIRECTORY, FILEINFO, etc */

#include "border.h"
    /* For: getBorderHeader(), etc */

#include "webcmiss.h"


#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>


/* -- Module Constants --------------------------------------------------------
*/

#define blort "blort"


/* -- Private Module Methods --------------------------------------------------
*/

void deleteAllFiles( CGI *cgi )
  {
  USER *user;
  JOB  *job;
	char *dir;
  FILEINFOITEM *item;
  FILEINFO     *file;
  char *fileName;
  char *dirName;

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
		writeError( "getJob() returned NULL in deleteAllFiles()\n" );
    disposeUser( user );
    return;
    }

	/* Now set about copying examples */
  dir = lookupString( cgi, "directory" );
  if( strncmp(dir,"input",5) == 0 )
    {
			item = job->inputFiles->fileInfoList->firstFileInfoItem;
			dirName = job->inputDirectoryName;
		}
	else
		{
			item = job->outputFiles->fileInfoList->firstFileInfoItem;
			dirName = job->outputDirectoryName;
		}

	while( item != NULL )
		{
    file = item->fileInfo;

		fileName = malloc( strlen( file->name ) + strlen( dirName ) + 2 );  /* + '/' + '\0' */

		if( fileName == NULL )
			{
			writeError( "Memory allocation failure in deleteAllFiles()\n" );
			return;
      }

    sprintf( fileName, "%s/%s", dirName, file->name );

		if (remove( fileName ))
			{
			writeError( "Error deleteing in deleteAllFiles()\n" );
			return;
			}

		free( fileName );
		item = item->nextFileInfoItem;
    }

	forwardUserDir( LIST_PROBLEM, user, dir );
  }



void printDeleteOKButton( FILE *output, USER *user, char *directory )
  {
  fprintf( output, "\n<FORM METHOD=\"POST\" ACTION=\"" DELETE_DIR "\">\n" );
	fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">\n",      user->name );
  fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">\n",  user->password );
  fprintf( output, "<INPUT NAME=\"directory\" TYPE=\"hidden\" VALUE=\"%s\">\n", directory );
  fprintf( output, "<INPUT NAME=\"stage\" TYPE=\"hidden\" VALUE=\"1\">\n" );
  fprintf( output, "<INPUT TYPE=\"submit\" VALUE=\"Delete Files\">\n" );
  fprintf( output, "</FORM>\n" );
  }


void printDeleteCancelButton( FILE *output, USER *user, char *directory )
  {
  fprintf( output, "\n<FORM METHOD=\"POST\" ACTION=\"" LIST_PROBLEM "\">\n" );
  fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->name );
  fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->password );
  fprintf( output, "<INPUT NAME=\"directory\" TYPE=\"hidden\" VALUE=\"%s\">\n", directory );
  fprintf( output, "<INPUT TYPE=\"submit\" VALUE=\"Cancel\">\n" );
  fprintf( output, "</FORM>\n" );
  }


void deleteFilesForm( CGI *cgi, FILE *output )
  {
  USER *user;
  JOB  *job;
	char *directory;

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
    writeError( "getJob() returned NULL in deleteFiles()\n" );
    disposeUser( user );
    return;
    }

  directory = lookupString( cgi, "directory" );
  if( directory == NULL )
    {    
    writeError( "No file to edit in deleteFiles()\n" );
    disposeUser( user );
    return;
    }

	fflush( stdout );  /* Hack to make sure stderr doesn't get seen first */

  fprintf( output, "<H1>Delete Files</H1>\n" );
  fprintf( output, "<P>" );
  fprintf( output, "<CENTER>" );
  fprintf( output, "<FONT SIZE=\"+1\">Are you sure you wish to delete the "
    "files in the %s directory?</FONT><P>\n", directory );
  fprintf( output, "</CENTER>" );

  fprintf( output, "<CENTER>" );
  fprintf( output, "<TABLE BORDER=\"0\" WIDTH=\"60%%\"><TR>" );
  fprintf( output, "<TD ALIGN=\"LEFT\">" );
  printDeleteOKButton( output, user, directory );
  fprintf( output, "</TD>" );
  fprintf( output, "<TD ALIGN=\"RIGHT\">" );
  printDeleteCancelButton( output, user, directory );
  fprintf( output, "</TD>" );
  fprintf( output, "</TR></TABLE>" );
  fprintf( output, "</CENTER>" );
  }



/* -- Program Entry Point -----------------------------------------------------
*/


int main( int argc, char *argv[] )
  {
	USER *user;
  CGI  *cgi;
	char *stagestring;
	int  stage;


  cgi = getCGIEnvironment( argc, argv );
  if( cgi == NULL )
    {
    printf( "<HR><P><B>ERROR:</B> getCGIEnvironment() failed in main()</P>"
            "<P>Email "WEBMASTER" <I>immediately</I></P>\n" );
    return 0;
    }

  user = getUser( cgi );
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
		/* Ask about deletion */
  	case 0:
			{
			borderGenerateHeader( stdout, user, False, NULL, NULL );
			deleteFilesForm( cgi, stdout );
			borderGenerateFooter( stdout );
			break;
	    }

		/* Do the deletion */
    case 1:
      {
			deleteAllFiles( cgi );
      break;
      }


    default:
      {
      writeError("Invalid stage number");
      break;
      }
    } 

  return 0;
  }
