/*
 *  Parameter-menu.c
 *
 *    Provides a menu for the user to select parameter examples.
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

#include "cgi-decode.h"
    /* For: CGI handling functions */

#include "border.h"
    /* For: getBorderHeader() */

#include "transfer.h"
    /* For: copy() */

#include "webcmiss.h"
    /* For: file location defs */



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


/* -- Parameter Menu ----------------------------------------------------------
*/

void printParameterMenu( FILE *output, USER *user )
  {
    fprintf( output, "<H1>CMISS Parameter Examples</H1>\n" );
    fprintf( output, "<P>\n" );
    fprintf( output, "<HR NOSHADE>\n" );
    fprintf( output, "<P>\n" );
    fprintf( output, "\n" );
    fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" PARAM_SET_EXAMPLE "\">\n" );
    fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->name );
    fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->password );

    catFile( output, PARAM_MENU_FILE );

    fprintf( output, "<P>\n" );
    fprintf( output, "<ol>\n" );
    fprintf( output, "<input type=\"submit\" value=\"Get Example\">\n" );
    fprintf( output, "</ol>\n" );
    fprintf( output, "</form></p>\n" );
  }


/* -- Program Entry Point -----------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
    USER *user;
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

    paramBorderGenerateHeader( stdout, user, False, NULL, NULL );
    printParameterMenu( stdout, user );
    borderGenerateFooter( stdout );

    return 0;
  }
