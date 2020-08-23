/*
 *  Border.c
 *
 *    Code for adding a border around the CGI generated web pages.
 */


/* -- Include Directives ------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stderr */

#include <stdlib.h>
    /* For: malloc(), free(), NULL, srand(), rand() */

#include <string.h>
    /* For: strcpy() */

#include "user-managment.h"
    /* For: USER* */

#include "border.h"
    /* For: our public methods and datatypes */

#include "webcmiss.h"
    /* For: statically defined file names etc */


static int borderCatFile( FILE *output, char *filename )
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

static void borderGenerateButton( FILE *output, USER *user, char *label, 
				  char *script, BOOLEAN showUser, char *dirLabel )
  {
    fprintf( output, "\t<TR>\n" );
    fprintf( output, "\t<TD CLASS=\"small\">\n" );
    fprintf( output, "\t<FORM METHOD=\"POST\" ACTION=\"%s\">\n",script );
    fprintf( output, "\t<INPUT TYPE=\"submit\" VALUE=\"%s\">\n",label );
    if ( showUser )
      {
	fprintf( output, "\t<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">\n",
		 user->name );
	fprintf( output, "\t<INPUT TYPE=\"hidden\" NAME=\"password\" VALUE=\"%s\">\n",
		 user->password );
      }
    if ( dirLabel )
      {
	fprintf( output, "\t<INPUT TYPE=\"hidden\" NAME=\"directory\" VALUE=\"%s\">\n",
		 dirLabel );
      }
    fprintf( output, "\t</FORM>\n" );
    fprintf( output, "\t</TD>\n" );
    fprintf( output, "\t</TR>\n" );
  }

static int borderGenerateSidebar( FILE *output, USER *user )
  {
    fprintf( output, "\t\n" );
    fprintf( output, "\t<TABLE CELLSPACING=\"0\" BORDER=\"0\" CELLPADDING=\"0\" WIDTH=\"100%%\">\n" );
    fprintf( output, "\t<TR VALIGN=\"top\">\n" );
    fprintf( output, "\t<TD BGCOLOR=\"#7a87ff\" WIDTH=\"10%%\">\n" );
    fprintf( output, "\t<TABLE CELLSPACING=\"0\" BORDER=\"0\" CELLPADDING=\"3\" WIDTH=\"100%%\">\n" );
    fprintf( output, "\t<TR>\n" );
    fprintf( output, "\t<TD>&nbsp;</TD>\n" );
    fprintf( output, "\t</TR>\n" );

    borderGenerateButton( output, user, "CMISS Login",     LOGIN_URL,       False,  NULL );
    borderGenerateButton( output, user, "Input Files",     LIST_PROBLEM,    True,   "input" );
    borderGenerateButton( output, user, "Output Files",    LIST_PROBLEM,    True,   "output" );
    borderGenerateButton( output, user, "Job Control",     JOB_CONTROL,     True,   NULL );
    borderGenerateButton( output, user, "Account Details", ACCOUNT_DETAILS, True,   NULL );

    fprintf( output, "\t<TR>\n" );
    fprintf( output, "\t<TD HEIGHT=\"650\">&nbsp;</TD>\n" );
    fprintf( output, "\t</TR>\n" );
    fprintf( output, "\t</TABLE>\n" );
    fprintf( output, "\t</TD>\n" );
    fprintf( output, "\t\n" );
    fprintf( output, "\t<TD VALIGN=\"top\">\n" );
    fprintf( output, "\t<TABLE CELLSPACING=\"0\" BORDER=\"0\" CELLPADDING=\"5\" WIDTH=\"100%%\">\n" );
    fprintf( output, "\t<TR>\n" );
    fprintf( output, "\t<TD VALIGN=\"top\">\n" );
    fprintf( output, "\t\n" );

    fflush( output );

    return 0;
  }

static int paramBorderGenerateSidebar( FILE *output, USER *user )
  {
    fprintf( output, "\t\n" );
    fprintf( output, "\t<TABLE CELLSPACING=\"0\" BORDER=\"0\" CELLPADDING=\"0\" WIDTH=\"100%%\">\n" );
    fprintf( output, "\t<TR VALIGN=\"top\">\n" );
    fprintf( output, "\t<TD BGCOLOR=\"#7a87ff\" WIDTH=\"10%%\">\n" );
    fprintf( output, "\t<TABLE CELLSPACING=\"0\" BORDER=\"0\" CELLPADDING=\"3\" WIDTH=\"100%%\">\n" );
    fprintf( output, "\t<TR>\n" );
    fprintf( output, "\t<TD>&nbsp;</TD>\n" );
    fprintf( output, "\t</TR>\n" );

    borderGenerateButton( output, user, "CMISS Login",     LOGIN_URL,   False,  NULL );
    borderGenerateButton( output, user, "Switch Examples", PARAM_MENU,  True,   "input" );

    fprintf( output, "\t<TR>\n" );
    fprintf( output, "\t<TD HEIGHT=\"650\">&nbsp;</TD>\n" );
    fprintf( output, "\t</TR>\n" );
    fprintf( output, "\t</TABLE>\n" );
    fprintf( output, "\t</TD>\n" );
    fprintf( output, "\t\n" );
    fprintf( output, "\t<TD VALIGN=\"top\">\n" );
    fprintf( output, "\t<TABLE CELLSPACING=\"0\" BORDER=\"0\" CELLPADDING=\"5\" WIDTH=\"100%%\">\n" );
    fprintf( output, "\t<TR>\n" );
    fprintf( output, "\t<TD VALIGN=\"top\">\n" );
    fprintf( output, "\t\n" );

    fflush( output );

    return 0;
  }

void writeHtmlHeader( FILE *output, char *titlestring, char *headingstring )
{
  fprintf( output, "<html><head>\n" );
  fprintf( output, "<title>%s</title>\n", titlestring );
  fprintf( output, "</head><body>\n" );
  fprintf( output, "<h1>%s</h1>\n", headingstring );
  fprintf( output, "<p></p><hr>\n" );
}



/* -- Public methods ---------------------------------------------------------
*/

void writeError(char *error_str)
{
  writeHtmlHeader( stdout,"CMISS error","Error:" );
  fprintf( stdout, "%s\n",error_str );
  fflush( stdout );
}

void borderWriteError( FILE *output, USER *user, char *error_str )
{
  borderGenerateHeader( output, user, False, NULL, NULL  );

  writeHtmlHeader( output, "CMISS error", "Error:" );
  fprintf( output, "%s\n", error_str );

  borderGenerateFooter( output );
}

int borderGenerateHeader( FILE *output, USER *user, BOOLEAN refresh,
                          char *address, char *other )
{
  int ierr;

  fprintf( output, "Content-type: text/html\n\n" );
  fprintf( output, "<HTML>\n" );
  fprintf( output, "<HEAD>\n" );
  if ( ierr = borderCatFile( output, BORDER_HEAD_FILE ) )
  {
    return ierr;
  }

  if ( refresh )
  {
    fprintf( output, "  <META HTTP-EQUIV=\"Refresh\" Content=\"10;"
      "URL=%s?user=%s&password=%s",
      address, user->name, user->password );
    if ( other != NULL )
	  {
	    fprintf( output, "&%s", other );
	  }
    fprintf( output, "\">\n" );
  }
  fprintf( output, "</HEAD>\n" );

  if ( ierr = borderCatFile( output, BORDER_HEADER_FILE ) )
  {
    return ierr;
  }

  return borderGenerateSidebar( output, user );
}

int borderGenerateFooter( FILE *output )
{
  return borderCatFile( output, BORDER_FOOTER_FILE );
}

void paramBorderWriteError( FILE *output, USER *user, char *error_str )
{
  paramBorderGenerateHeader( output, user, False, NULL, NULL );

  writeHtmlHeader( output, "CMISS error", "Error:" );
  fprintf( output, "%s\n", error_str );

  borderGenerateFooter( output );
}

int paramBorderGenerateHeader( FILE *output, USER *user, BOOLEAN refresh,
                               char *address, char *other )
{
  int ierr;

  fprintf( output, "Content-type: text/html\n\n" );
  fprintf( output, "<HTML>\n" );
  fprintf( output, "<HEAD>\n" );
  if ( ierr = borderCatFile( output, BORDER_HEAD_FILE ) )
  {
    return ierr;
  }

  if ( refresh )
  {
    fprintf( output, "  <META HTTP-EQUIV=\"Refresh\" Content=\"10;"
      "URL=%s?user=%s&password=%s",
      address, user->name, user->password );
    if ( other != NULL )
	  {
	    fprintf( output, "&%s", other );
	  }
    fprintf( output, "\">\n" );
  }
  fprintf( output, "</HEAD>\n" );

  if ( ierr = borderCatFile( output, BORDER_HEADER_FILE ) )
  {
    return ierr;
  }

  return paramBorderGenerateSidebar( output, user );
}

