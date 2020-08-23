/*
 *  Delete File.c
 *
 *    Deletes a file from the Output directory.
 *
 *    "We come from a blue planet light years away,
 *     where everything multiplies at an amazing rate,
 *     we're out here in the universe buying real estate,
 *     hope we havn't gotten here to late.
 *
 *     We're humans from Earth.
 *     We're humans from Earth.
 *     You have nothing at all to fear,
 *     I think we're going to like it here.
 *
 *     We're looking for a planet with atmosphere,
 *     where the air is fresh and the water clear,
 *     with lots of sun like you have here,
 *     three or four hundred days a year.
 *
 *     We're humans from Earth.
 *
 *     Bought manhatten for a string of beads. 
 *     We bought along some gadgets for you to see.
 *     Here a crazy little thing we call TV.
 *     Do you have electricity?
 *
 *     We're humans from Earth.
 *     
 *     I know we may seem pretty strange to you,
 *     but we got know-how and golden rules.
 *     We're here to see manifest destiny through.
 *     There ain't nothing we can't get used to.
 *
 *     We're humans from Earth."
 *
 *                                     -- Humans From Earth, T-Bone Burnett
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


/* -- Module Constants --------------------------------------------------------
*/

#define blort "blort"


/* -- Deletes ----------------------------------------------------------
*/

void deleteFile( CGI *cgi )
  {
  USER *user;
  JOB  *job;
  char *fileName;
  char *fullFileName;
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
		writeError( "getJob() returned NULL in deleteFile()\n" );
    disposeUser( user );
    return;
    }
  
	/* Now set about copying examples */
  directory = lookupString( cgi, "directory" );
  if( directory == NULL )
    {
    directory = "input";
		}

  fileName = lookupString( cgi, "file" );
  if( fileName == NULL )
    {
		writeError( "No file name specified in deleteFile()\n" );
    return;
    }
  
  if( strchr( fileName, '/' ) != NULL || strstr( "..", fileName ) != NULL )
    {
		writeError( "Illegal character in filename in deleteFile()\n" );
    return; 
    }

  fullFileName = malloc( strlen( fileName ) + strlen( job->outputDirectoryName )
    + 2 );  /* + '/' + '\0' */
  if( fullFileName == NULL )
    {
		writeError( "Memory allocation failure in deleteFile()\n" );
    return;
    }

  sprintf( fullFileName, "%s/%s", job->outputDirectoryName, fileName );

  remove( fullFileName );

  forwardUserDir( LIST_PROBLEM, user, directory );
  }


/* -- Asks confirmation to delete --------------------------------------------------
*/

void printDeleteOKButton( FILE *output, USER *user, char *file, char *directory )
  {
  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" DELETE_FILE "\">\n" );
  fprintf( output, "<INPUT NAME=\"file\" TYPE=\"hidden\" VALUE=\"%s\">\n", file );
  fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->name );
  fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->password );
  fprintf( output, "<INPUT NAME=\"directory\" TYPE=\"hidden\" VALUE=\"%s\">\n", directory );
  fprintf( output, "<INPUT NAME=\"stage\" TYPE=\"hidden\" VALUE=\"1\">\n" );
  fprintf( output, "<INPUT TYPE=\"submit\" VALUE=\"Delete File\">\n" );
  fprintf( output, "</FORM>\n" );
  }


void printDeleteCancelButton( FILE *output, USER *user, char *directory )
  {
  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" LIST_PROBLEM "\">\n" );
  fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->name );
  fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->password );
  fprintf( output, "<INPUT TYPE=\"submit\" VALUE=\"Cancel\">\n" );
  fprintf( output, "<INPUT NAME=\"directory\" TYPE=\"hidden\" VALUE=\"%s\">\n", directory );
  fprintf( output, "</FORM>\n" );
  }


void deleteFileForm( CGI *cgi, FILE *output )
  {
  USER *user;
  JOB  *job;
  char *fileName;
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
    fprintf( output, "getJob() returned NULL in revertFile()\n" );
    disposeUser( user );
    return;
    }

  fileName = lookupString( cgi, "file" );
  if( fileName == NULL )
    {    
    fprintf( output, "No file to delete in deleteFile()\n" );
    disposeUser( user );
    return;
    } 

  if( strchr( fileName, '/' ) != NULL || strstr( "..", fileName ) != NULL )
    {
    fprintf( output, "Illegal character in filename in revertFile()\n" );
    disposeUser( user );
    return; 
    }

  directory = lookupString( cgi, "directory" );
  if( directory == NULL )
    {    
    directory = "output";
    }

  fprintf( output, "<HR NOSHADE><P>\n" );
  printJobHeading( output, user, job );
  fprintf( output, "<P>\n" );

  fprintf( output, "<CENTER>\n" );
  fprintf( output, "<FONT SIZE=\"+1\">Are you sure you wish to delete the "
    "file <B>%s</B>?</FONT><P>\n", fileName );
  fprintf( output, "</CENTER>\n" );

  fprintf( output, "<CENTER>\n" );
  fprintf( output, "<TABLE BORDER=\"0\" WIDTH=\"60%%\"><TR>\n" );
  fprintf( output, "<TD ALIGN=\"LEFT\">\n" );
  printDeleteOKButton( output, user, fileName, directory );
  fprintf( output, "</TD>\n" );
  fprintf( output, "<TD ALIGN=\"RIGHT\">\n" );
  printDeleteCancelButton( output, user, directory );
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
			deleteFileForm( cgi, stdout );
			borderGenerateFooter( stdout );
			break;
	    }

		/* Do the deletion */
    case 1:
      {
			deleteFile( cgi );
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
