/*
 *  Revert File.c
 *
 *    Revert File dialog for the CMISS web interface
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stdout, FILE, etc */

#include <stdlib.h>
    /* For: malloc(), free(), NULL */

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


/* -- Private Module Methods --------------------------------------------------
*/



/* -- Revert the file -------------------------------------------------
*/

BOOLEAN fileOK( char *fileName )
  {
  FILEINFO *fileInfo;

  fileInfo = newFileInfo( fileName );
  if( fileInfo == NULL )
    {
    return False;
    }

  disposeFileInfo( fileInfo );

  return True;
  }


void copyFile( char *oldFile, char *newFile )
  {
  FILE *source;
  FILE *destination;
  char  buffer[512];

  source = fopen( oldFile, "rt" );
  if( source == NULL )
    {
/* Error message here later */
    return;
    }

  destination = fopen( newFile, "wt" );
  if( destination == NULL )
    {
/* Error message here later */
    return;
    }

  fgets( buffer, 512, source );
  while( !feof( source ) )
    {
    fputs( buffer, destination );
    fgets( buffer, 512, source );
    }
  }



/* -- Module Methods ----------------------------------------------------------
*/

void revertFile( CGI *cgi )
  {
  USER *user;
  JOB  *job;
  char *fileName;
  char *newFile;
  char *oldFile;
#ifdef SETGID
	struct group *cmiss;
#endif

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
    fprintf( stderr, "getJob() returned NULL in revertFile()\n" );
    disposeUser( user );
    return;
    }
  
#ifdef SETGID
	/* Set GID to cmiss */
	if ((cmiss = getgrnam("cmiss")) == NULL)
		{
    fprintf( stderr, "getgrnam() error in revertFile()\n" );
    disposeUser( user );
    return;
		}
	if (setgid(cmiss->gr_gid))
		{
    fprintf( stderr, "setgid() error in revertFile()\n" );
    disposeUser( user );
    return;
		}
#endif

  fileName = lookupString( cgi, "file" );
  if( fileName == NULL )
    {
		writeError( "No file name specified in revertFile()\n" );
    return;
    }
  
  if( strchr( fileName, '/' ) != NULL || strstr( "..", fileName ) != NULL )
    {
		writeError( "Illegal character in filename in revertFile()\n" );
    return; 
    }

  newFile = malloc( strlen( fileName ) + strlen( job->inputDirectoryName )
    + 2 );  /* + '/' + '\0' */
  if( newFile == NULL )
    {
		writeError( "Memory allocation failure in revertFile() [1]\n" );
    return;
    }
  
  sprintf( newFile, "%s/%s", job->inputDirectoryName, fileName );

  oldFile = malloc( strlen( fileName ) + strlen( job->originalDirectoryName )
    + 2 );  /* + '/' + '\0' */
  if( oldFile == NULL )
    {
		writeError( "Memory allocation failure in revertFile() [2]\n" );
    return;
    }

  sprintf( oldFile, "%s/%s", job->originalDirectoryName, fileName );

  if( fileOK( newFile ) == False || fileOK( oldFile ) == False )
    {
		writeError( "A problem occured with the files.\n" );
    return;
    }

  copyFile( oldFile, newFile );

	forwardUser( LIST_PROBLEM, user );
  }


/* -- Print the revert file form -------------------------------------------------
*/


void printRevertOKButton( FILE *output, USER *user, char *file )
  {
  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" REVERT_FILECONF "\">\n" );
  fprintf( output, "<INPUT NAME=\"file\" TYPE=\"hidden\" VALUE=\"%s\">", file );
  fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">", user->name );
  fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">", user->password );
  fprintf( output, "<INPUT NAME=\"stage\" TYPE=\"hidden\" VALUE=\"1\">\n" );
  fprintf( output, "<INPUT TYPE=\"submit\" VALUE=\"Revert File\">\n" );
  }


void printRevertCancelButton( FILE *output, USER *user )
  {
  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" LIST_PROBLEM "\">\n" );
  fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->name );
  fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->password );
  fprintf( output, "<INPUT TYPE=\"submit\" VALUE=\"Cancel\">\n" );
  fprintf( output, "</FORM>\n" );
  }


void revertFileForm( CGI *cgi, FILE *output )
  {
  USER *user;
  JOB  *job;
  char *fileName;
  char *fullChangedName;
  char *fullOriginalName;

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
    fprintf( stderr, "getJob() returned NULL in revertFile()\n" );
    disposeUser( user );
    return;
    }

  fileName = lookupString( cgi, "file" );
  if( fileName == NULL )
    {    
    fprintf( stderr, "No file to edit in revertFile()\n" );
    disposeUser( user );
    return;
    } 

  if( strchr( fileName, '/' ) != NULL || strstr( "..", fileName ) != NULL )
    {
    fprintf( stderr, "Illegal character in filename in revertFile()\n" );
    disposeUser( user );
    return; 
    }

  fullChangedName = malloc( strlen( job->inputDirectoryName ) +
    strlen( fileName ) + 2 ); /* + '/' + '\0' */
  if( fullChangedName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in listProblem()\n" );
    disposeUser( user );
    return;
    }
  sprintf( fullChangedName, "%s/%s", job->inputDirectoryName, fileName );

  fullOriginalName = malloc( strlen( job->originalDirectoryName ) + 
    strlen( fileName ) + 2 );  /* + '/' + '\0' */
  if( fullOriginalName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in listProblem()\n" );
    disposeUser( user );
    return;
    }
  sprintf( fullOriginalName, "%s/%s", job->originalDirectoryName,
    fileName );

  fprintf( output, "<HR NOSHADE><P>\n" );
  printJobHeading( output, user, job );
  fprintf( output, "<P>\n" );

  fprintf( output, "<FONT SIZE=\"+1\">\n" );
	fprintf( output, "Are you sure you wish to revert the file <B>%s</B>\n", fileName );
	fprintf( output, "back to the original example version (discarding all changes\n" );
	fprintf( output, "you have made to it)?\n" );
	fprintf( output, "</FONT><P>\n" );

  fprintf( output, "<CENTER>\n" );

  fprintf( output, "<TABLE BORDER=\"0\" WIDTH=\"60%%\"><TR>\n" );
  fprintf( output, "<TD ALIGN=\"LEFT\">\n" );
  printRevertOKButton( output, user, fileName );
  fprintf( output, "</TD>\n" );
  fprintf( output, "<TD ALIGN=\"RIGHT\">\n" );
  printRevertCancelButton( output, user );
  fprintf( output, "</TD>\n" );
  fprintf( output, "</TR></TABLE>\n" );
  fprintf( output, "</CENTER>\n" );
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
            "<P>Email "WEBMASTER"<I>immediately</I></P>\n" );
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
			revertFileForm( cgi, stdout );
			borderGenerateFooter( stdout );
			break;
	    }

		/* Do the deletion */
    case 1:
      {
			revertFile( cgi );
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
