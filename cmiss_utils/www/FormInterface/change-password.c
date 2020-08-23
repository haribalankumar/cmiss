/*
 *  Change Password.c
 *
 *    Changes the users password.
 */


/* -- Include Files -----------------------------------------------------------
*/


#include <stdio.h>
    /* For: fprintf(), stderr, etc */

#include <stdlib.h>
    /* For: malloc(), free(), NULL, etc */

#include <string.h>
    /* For: strlen(), strcpy(), strstr(), etc */

#include <crypt.h>
    /* For: crypt() */

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
    /* For: borderGenerateHeader() */

#include "webcmiss.h"

#include "password.h"


/* -- Module Datatypes --------------------------------------------------------
*/

#if !defined( BOOLEAN )
#define BOOLEAN int
#define True 1
#define False 0
#endif

/* -- Module Constants --------------------------------------------------------
*/



/* -- Module Private Methods --------------------------------------------------
*/ 


/* -- Module Methods ----------------------------------------------------------
*/

void changePasswd( FILE *output, CGI *cgi )
  {
  USER      *user;
  JOB       *job;
	char      *oldpasswd;
	char      *newpasswd1;
	char      *newpasswd2;
	char      *password;

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
    fprintf( stderr, "getJob() returned NULL in changePasswd()\n" );
    disposeUser( user );
    return;
    }

	/* If the password change has been canceled, then return */
	if( lookupString( cgi, "cancel" ))
		{
		forwardUser( ACCOUNT_DETAILS, user ); 
    disposeUser( user );
    disposeJob( job );
    return;
		}

	/* Get the oldpassword, and the two copies of the new one */
	oldpasswd  = lookupString( cgi, "oldpassword" );
	newpasswd1 = lookupString( cgi, "newpassword1" );
	newpasswd2 = lookupString( cgi, "newpassword2" );
  if( oldpasswd == NULL || newpasswd1 == NULL || newpasswd2 == NULL )
    {
		borderGenerateHeader( output, user, False, NULL, NULL );
    fprintf( output, "All password forms not filled in." );
    disposeUser( user );
    disposeJob( job );
		borderGenerateFooter( output );
    return;
    }

	/* Check the old password is valid */
	if (strcmp( crypt( oldpasswd, user->name ), user->password ) != 0)	{
		borderGenerateHeader( output, user, False, NULL, NULL );
    fprintf( output, "Username or password incorrect." );
    disposeUser( user );
    disposeJob( job );
		borderGenerateFooter( output );
    return;
	}

	/* Check length of the password */
	if (strlen(newpasswd1) < 6) {
		borderGenerateHeader( output, user, False, NULL, NULL );
    fprintf( output, "New password should be at least 6 characters long." );
    disposeUser( user );
    disposeJob( job );
		borderGenerateFooter( output );
    return;
	}

	/* Check the two new passwords are the same */
	if (strcmp(newpasswd1,newpasswd2) != 0)	{
		borderGenerateHeader( output, user, False, NULL, NULL );
    fprintf( output, "New passwords do not match." );
    disposeUser( user );
    disposeJob( job );
		borderGenerateFooter( output );
    return;
	}

	/* Encrypt it */
	password = crypt( newpasswd1, user->name );

	/* Change the password */
	if (changePassword( user->name, password )) {
		borderGenerateHeader( output, user, False, NULL, NULL );
    fprintf( output, "Error setting the new password." );
    disposeUser( user );
    disposeJob( job );
		borderGenerateFooter( output );
    return;
	}
  user->password = password;
	/* 
		 strdup(password);
		 stripTrailingSpace( user->password );
		 */

	/* Onto the file lister .. */
  forwardUser( ACCOUNT_DETAILS, user ); 
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
            "<P>Email "WEBMASTER"<I>immediately</I></P>\n" );
    return 0;
    }

  changePasswd( stdout, cgi );

  return 0;
  }


