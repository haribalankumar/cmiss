/*
 *  Check Pasword.c
 *
 *    Perform the initial check of the user password, and encrypt it 
 *    for further use. The user is then passed onto list-problem.cgi,
 *    or parameter-menu.cgi.
 *
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), FILE */

#include "user-managment.h"
    /* For: USER, USERLIST, readUserList(), passwordCheck(), forwardUser(), etc */

#include "get-user.h"
    /* For: getUser() */

#include "webcmiss.h"


/* -- Module Methods ----------------------------------------------------------
*/

void checkPassword( CGI *cgi )
  {
  USER *user;
  char *system;

  user = getUserCrypt( cgi );
  if( user == NULL )
    {
    /* Not always an error, the user could simply have supplied an incorrect
       password. getUser() will report back those conditions to the user
       (one hopes...) */
    return;
    }

  /* Check to see if we have a variable "system" defined as "parameter". If
     so use the simplified parameter system, otherwise use the default. */
  system = lookupString( cgi, "system" );
  if( system != NULL )
    {
    if( strcmp(system,"parameter") == 0 )
      {
	forwardUser( PARAM_MENU, user );
	return;
      }
    }

   forwardUser( LIST_PROBLEM, user );
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

  checkPassword( cgi );

  return 0;
  }
