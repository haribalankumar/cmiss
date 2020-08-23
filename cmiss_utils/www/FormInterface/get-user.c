/*
 *  Get User.c
 *
 *    Checks and gets a user record given a CGI structure
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stderr */

#include <stdlib.h>
    /* For: NULL, malloc(), free() */

#include <crypt.h>
    /* For: crypt() */

#include "user-managment.h"
    /* For: USER, USERLIST, readUserList(), lookupUser(), passwordCheck() */

#include "cgi-decode.h"
    /* For: CGI handling functions */

#include "get-user.h"
    /* For: our public interface */

#include "tags.h"
    /* For: tagCompare() */

#include "print-error.h"
    /* For: printErrorMessage() */

#include "webcmiss.h"


/* Password Constants */
typedef enum
  {
  PasswordInvalid = 0,
  PasswordValid   = 1
  }
PasswordStatus;

#if !defined( BOOLEAN )
#define BOOLEAN int
#define True 1
#define False 0
#endif



/* -- Module Private Methods --------------------------------------------------
*/

void displayNoPasswordDialog( void )
  {
  fprintf( stdout, "<P>You have no username or password<P>" );

  fprintf( stdout, "<FORM METHOD=\"POST\" ACTION=\"password.cgi\">" );
  }


void displayIncorrectPasswordDialog( void )
  {
    fprintf( stdout, "<P>The username or password you entered is incorrect<P>" );
  }


PasswordStatus passwordCheck( USER *user, char *password, BOOLEAN first )
  {
	if (first)
		{
    if( tagCompare( user->password, crypt( password, user->name ) ) == True )
      return PasswordValid;
    else
      return PasswordInvalid;
		}
	else
		{
		if( tagCompare( user->password, password ) == True )
			return PasswordValid;
		else
			return PasswordInvalid;
		}
  }


USER *getUserPrivate( CGI *cgi, BOOLEAN first )
  {
  USERLIST  *userList;
  USER      *user;
  USER      *userCopy;
  char      *name;
  char      *password;

  userList = newUserList();
  if( userList == NULL )
    {
    printErrorMessage( "Login Error", "newUserList() returned NULL in getUser()" );
    /* fprintf( stderr, "newUserList() returned NULL in getUser()\n" ); */
    return NULL;
    }

  readUserList( userList, ESU_MASTER_USER_LIST );
  
  name = lookupString( cgi, "user" );
  password = lookupString( cgi, "password" );

  if( name == NULL )
    {
    printErrorMessage( "Login Error", "Username or Password not provided" );
    /* displayNoPasswordDialog(); */
    disposeUserList( userList );
    return NULL;
    }

  user = lookupUser( name, userList );
  if( user == NULL )
    {
	  printErrorMessage( "Login Error", "Username or Password invalid" );
    /* displayIncorrectPasswordDialog(); */
    disposeUserList( userList );
    return NULL;
    }
  
  if( passwordCheck( user, password, first ) == PasswordInvalid )
    {
	  printErrorMessage( "Login Error", "Username or Password invalid" );
    /* displayIncorrectPasswordDialog(); */
    disposeUserList( userList );
    return NULL;
    }

  userCopy = duplicateUser( user );
  if( userCopy == NULL )
    {
    fprintf( stderr, "duplicateUser() returned NULL in getUser()\n" );
    disposeUserList( userList );
    return NULL;
    }

  disposeUserList( userList );

  return userCopy;
  }


/* -- Public Methods ----------------------------------------------------------
*/

USER *getUser( CGI *cgi )
  {
  return getUserPrivate( cgi, False );
  }


USER *getUserCrypt( CGI *cgi )
  {
  return getUserPrivate( cgi, True );
  }

