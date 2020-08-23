/*
 *  Parameter Example.c
 *
 *    Displays the parameter example contained in the users directory.
 */

/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stderr, etc */

#include <stdlib.h>
    /* For: malloc(), free(), NULL, etc */

#include "get-user.h"
    /* For: getUser() */

#include "user-managment.h"
    /* For: USER, disposeUser() */

#include "job.h"
    /* For: JOB */

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
    strcat( fullName, "/" );
    strcat( fullName, subname );

    return fullName;
  }


static int catFile( FILE *output, char *filename )
  {
    FILE *input;
    int c;

    if ( ( input = fopen( filename, "r" ) ) == NULL )
      {
	fprintf( output, "borderCatFile: Cannot open file \"%s\"\n", filename );
	return -1;
      }

    while ((c = fgetc(input)) != EOF)
      {
	fputc(c, output);
      }

    fclose(input);
    fflush( output );

    return 0;
  }


void printSubmitButton( FILE *output, USER *user, char *example )
  {
  /* fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" PARAM_JOB "\">" ); */
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">",  user->name );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"password\" VALUE=\"%s\">", user->password );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"example\" VALUE=\"%s\">", example );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"stage\" VALUE=\"0\">" );
  fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"Submit Job\">" );
  /* fprintf( output, "</FORM>" ); */
  }


void printResetButton( FILE *output, USER *user, char *example )
  {
  /* fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" PARAM_JOB "\">" ); */
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">",  user->name );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"password\" VALUE=\"%s\">", user->password );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"example\" VALUE=\"%s\">", example );
  fprintf( output, "<INPUT TYPE=\"reset\" VALUE=\"Reset Values\">" );
  /* fprintf( output, "</FORM>" ); */
  }


void printAbortButton( FILE *output, USER *user, char *example )
  {
  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" PARAM_JOB "\">" );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">",  user->name );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"password\" VALUE=\"%s\">", user->password );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"example\" VALUE=\"%s\">", example );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"stage\" VALUE=\"10\">" );
  fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"Abort Job\">" );
  fprintf( output, "</FORM>" );
  }


void printViewButton( FILE *output, USER *user, char *example )
  {
  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" PARAM_JOB "\">" );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">",  user->name );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"password\" VALUE=\"%s\">", user->password );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"example\" VALUE=\"%s\">", example );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"stage\" VALUE=\"20\">" );
  fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"View Results\">" );
  fprintf( output, "</FORM>" );
  }


void printParameterExample( FILE *output, USER *user, JOB *job, char *example )
  {
    catFile( output, generateSubName( user->homeDirectory, "Input/param.html" ) );

    fprintf( output, "<P>\n" );
    printBriefJobHeading( output, job );
    fprintf( output, "<P>\n" );

    fprintf( output, "<CENTER>\n" );
    fprintf( output, "<TABLE BORDER=\"0\" WIDTH=\"30%%\">\n" );
    fprintf( output, "<TR>\n" );
    fprintf( output, "<TD ALIGN=\"LEFT\" WIDTH=\"20%%\">\n" );
    printSubmitButton( output, user, example );
    fprintf( output, "</TD>\n" );
    fprintf( output, "<TD ALIGN=\"LEFT\" WIDTH=\"20%%\">\n" );
    printResetButton( output, user, example );
    fprintf( output, "</TD>\n" );
    fprintf( output, "</TR>\n" );
    fprintf( output, "</TABLE>\n" );
    fprintf( output, "</CENTER>\n" );
    fprintf( output, "</FORM>" );
    fprintf( output, "<P>\n" );

    fprintf( output, "<CENTER>\n" );
    fprintf( output, "<TABLE BORDER=\"0\" WIDTH=\"30%%\">\n" );
    fprintf( output, "<TR>\n" );
    fprintf( output, "<TD ALIGN=\"LEFT\" WIDTH=\"20%%\">\n" );
    printAbortButton( output, user, example );
    fprintf( output, "</TD>\n" );
    fprintf( output, "<TD ALIGN=\"LEFT\" WIDTH=\"20%%\">\n" );
    printViewButton( output, user, example );
    fprintf( output, "</TD>\n" );
    fprintf( output, "</TR>\n" );
    fprintf( output, "</TABLE>\n" );
    fprintf( output, "</CENTER>\n" );

    fprintf( output, "<P>\n" );
  }


/* -- Program Entry Point -----------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
    USER *user;
    CGI   *cgi;
    JOB   *job;
    char *buff;
    char *example;


    cgi = getCGIEnvironment( argc, argv );
    if( cgi == NULL ) {
      printf( "<HR><P><B>ERROR:</B> getCGIEnvironment() failed in main()</P>"
        "<P>Email "WEBMASTER"<I>immediately</I></P>\n" );
      return -1;
    }

    user = getUser( cgi );
    if( user == NULL ) {
      /* Not always an error, the user could simply have supplied an incorrect
         password. getUser() will report back those conditions to the user
         (one hopes...) */
      return -1;
    }

    job = getJob( user );
    if( job == NULL ) {
      paramBorderWriteError( stdout, user, "No job in paramExample()\n" );
      disposeUser( user );
      return -1;
    }

    /* Now set about copying examples */
    example = lookupString( cgi, "example" );
    if( example == NULL ) {
      paramBorderWriteError( stdout, user, "No example name in paramExample()\n" );
      return -1;
    }

    if ( strncmp( job->status, "running", 7 ) == 0 ) {
      buff = malloc( strlen( example ) + 16 );
      if ( buff == NULL ) {
        paramBorderWriteError( stdout, user, "Malloc failure in paramExample()\n" );
        return -1;
      }
      
      sprintf( buff, "example=%s", example );
      paramBorderGenerateHeader( stdout, user, True, PARAM_EXAMPLE, buff );
      free( buff );
    }
    else {
      paramBorderGenerateHeader( stdout, user, False, NULL, NULL );
    }
    printParameterExample( stdout, user, job, example );
    borderGenerateFooter( stdout );

    return 0;
  }
