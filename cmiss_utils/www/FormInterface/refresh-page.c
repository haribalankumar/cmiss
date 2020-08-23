/*
 *  Refresh Page.c
 *
 *    Refreshes a page. By using this cgi program we avoid those ugly 
 *    get URLs.
 *
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), FILE */

#include <stdlib.h>
    /* For: malloc, free, putenv */

#include "user-managment.h"
    /* For: USER, USERLIST, readUserList(), passwordCheck(), forwardUser(), etc */

#include "get-user.h"
    /* For: getUser() */

#include "webcmiss.h"


/* -- Module Methods ----------------------------------------------------------
*/

void refreshPage( CGI *cgi )
  {
  USER *user;
  char *stage;
  char *address;
  char *buff;

  user = getUser( cgi );
  if( user == NULL )
    {
    /* Not always an error, the user could simply have supplied an incorrect
       password. getUser() will report back those conditions to the user
       (one hopes...) */
    return;
    }

  /* Find the address of the refering page */
  address = lookupString( cgi, "address" );
  if( address == NULL )
    {
      borderWriteError( stdout, user, "Error finding address" );
      return;
    }

  /* Check to see if we have a variable "stage" -- if so this should be set. */
  stage = lookupString( cgi, "stage" );
  if( stage != NULL )
    {
      buff = malloc( strlen( stage ) + 8 );
      if ( buff == NULL )
	{
	  return;
	}

      sprintf( buff, "stage=%s", stage );
      forwardUserTag( address, user, buff );
      free( buff );
      return;
    }

  forwardUser( address, user );
  }



/* -- Program Entry Point -----------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
  CGI *cgi;

  cgi = getCGIEnvironment( argc, argv );
  if( cgi == NULL )
    {
    printf( "<HR><P><B>ERROR:</B> getCGIEnvironment() failed in main()</P>"
            "<P>Email "WEBMASTER" <I>immediately</I></P>\n" );
    return 0;
    }

  refreshPage( cgi );

  return 0;
  }
