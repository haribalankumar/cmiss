/*
 *  User Managment.c
 *
 *    Manages Users
 */


/* -- Include Directives ------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stderr */

#include <stdlib.h>
    /* For: malloc(), free(), NULL, srand(), rand() */

#include <string.h>
    /* For: strcpy() */

#include <ctype.h>
    /* For: isprint(), tolower() */

#include <time.h>
    /* For: time() */

#include <unistd.h>
    /* For: exec() */

#include "groups.h"
    /* For: GROUP, GROUPFILE, setGroupFilePos(), etc, etc */

#include "rules.h"
    /* For: RULES, newRules(), etc */

#include "rule-parse.h"
    /* For: parseRules() */

#include "print-error.h"
    /* For: parseErrorMessage() */

#include "user-managment.h"
    /* For: our public methods and datatypes */



/* -- Private Module Constants ------------------------------------------------
*/

/* String Constants */
#define StringMatch 0



/* -- Private Module Utilities ------------------------------------------------
*/

void stripTrailingSpace( char *string )
  {
  char *p;
  char *lastSpace;

  p = string + strlen( string );
  lastSpace = p;

  while( p != string )
    {
    if( *p == '\r' || *p == '\n' || *p == '\t' || *p == ' ' || *p == '\0' )
      lastSpace = p;
    else
      break;

    p--;
    }

  *lastSpace = '\0';
  }



/* -- Public Ctors and Dtors --------------------------------------------------
*/

USER *newUser( char *name, char *fullName, char *password, char *homeDirectory,
  char *email )
  {
  USER *tempUser;

  tempUser = malloc( sizeof( USER ) );
  if( tempUser == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newUser() [1]\n" );
    return NULL;
    }

  tempUser->name = malloc( strlen( name ) + 1 );
  if( tempUser->name == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newUser() [2]\n" );
    free( tempUser );
    return NULL;
    }

  tempUser->fullName = malloc( strlen( fullName ) + 1 );
  if( tempUser->fullName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newUser() [3]\n" );
    free( tempUser->name );
    free( tempUser );
    return NULL;
    }

  tempUser->password = malloc( strlen( password ) + 1 );
  if( tempUser->password == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newUser() [4]\n" );
    free( tempUser->fullName );
    free( tempUser->name );
    free( tempUser );
    return NULL;
    }

  tempUser->homeDirectory = malloc( strlen( homeDirectory ) + 1 );
  if( tempUser->homeDirectory == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newUser() [5]\n" );
    free( tempUser->password );
    free( tempUser->fullName );
    free( tempUser->name );
    free( tempUser );
    return NULL;
    }

  tempUser->email = malloc( strlen( email ) + 1 );
  if( tempUser->email == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newUser() [6]\n" );
    free( tempUser->homeDirectory );
    free( tempUser->password );
    free( tempUser->fullName );
    free( tempUser->name );
    free( tempUser );
    return NULL;
    }

  strcpy( tempUser->name, name );
  strcpy( tempUser->fullName, fullName );
  strcpy( tempUser->password, password );
  strcpy( tempUser->homeDirectory, homeDirectory );
  strcpy( tempUser->email, email );

  stripTrailingSpace( tempUser->name );
  stripTrailingSpace( tempUser->fullName );
  stripTrailingSpace( tempUser->password );
  stripTrailingSpace( tempUser->homeDirectory );
  stripTrailingSpace( tempUser->email );

  return tempUser;
  }


USER *duplicateUser( USER *user )
  {
  USER *tempUser;

  tempUser = newUser( user->name, user->fullName, user->password, 
										  user->homeDirectory, user->email );
  if( tempUser == NULL )
    {
    fprintf( stdout, "newUser() returned NULL in duplicateUser()\n" );
    return NULL;
    }

  return tempUser;
  }


void disposeUser( USER *user )
  {
  if( user->name != NULL )
    free( user->name );
  if( user->fullName != NULL )
    free( user->fullName );
  if( user->password != NULL )
    free( user->password );
  if( user->homeDirectory != NULL )
    free( user->homeDirectory );
  if( user->email != NULL )
    free( user->email );
  free( user );
  }


USERITEM *newUserItem( USER *user )
  {
  USERITEM *tempUserItem;
  
  tempUserItem = malloc( sizeof( USERITEM ) );
  if( tempUserItem == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newUserItem()\n" );
    return NULL;
    }

  tempUserItem->user = user;
  tempUserItem->nextUserItem = NULL;

  return tempUserItem;
  }


void disposeUserItem( USERITEM *userItem )
  {
  disposeUser( userItem->user );
  free( userItem );
  }


void disposeUserItems( USERITEM *userItem )
  {
  if( userItem != NULL )
    {
    disposeUserItems( userItem->nextUserItem );
    disposeUserItem( userItem );
    }
  }


USERLIST *newUserList( void )
  {
  USERLIST *tempUserList;

  tempUserList = malloc( sizeof( USERLIST ) );
  if( tempUserList == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newUserList()\n" );
    return NULL;
    }

  tempUserList->firstUserItem = NULL;

  return tempUserList;
  }


void disposeUserList( USERLIST *userList )
  {
  if( userList->firstUserItem != NULL )
    disposeUserItems( userList->firstUserItem );

  free( userList );
  }



/* -- User List Input ---------------------------------------------------------
*/

USER *lookupUser( char *name, USERLIST *userList )
  {
  USERITEM *userItem;
  USER     *user;

  if( name == NULL )
    return NULL;

  userItem = userList->firstUserItem;
  while( userItem != NULL )
    {
    user = userItem->user;

    if( tagCompare( name, user->name ) == True )
      return user;
    
    userItem = userItem->nextUserItem;
    }

  return NULL;
  }


void addUserToUserList( USER *user, USERLIST *userList )
  {
  USERITEM *userItem, *userPtr;

  userItem = newUserItem( user );
  if( userItem == NULL )
    {
    fprintf( stderr, "newUserItem() returned NULL in addUserToUserList()\n" );
    return;
    }

#ifdef OLD_CODE
	/* This code reverses the file order */
  userItem->nextUserItem = userList->firstUserItem;
  userList->firstUserItem = userItem;
#else
	/* This code doesn't reverse the file order */
	userItem->nextUserItem = NULL;

	if (userList->firstUserItem == NULL) {
		userList->firstUserItem = userItem;
	}
	else {
		userPtr = userList->firstUserItem;
		while( userPtr->nextUserItem != NULL )
			{
			userPtr = userPtr->nextUserItem;
			}
		userPtr->nextUserItem = userItem;
	}
#endif
  }



void readUserList( USERLIST *userList, char *baseName )
  {
  GROUPFILE *file;
  RULES     *rules;
  GROUP     *group;
  TAG       *userTag;
  TAG       *nameTag;
  TAG       *passwordTag;
  TAG       *homeDirectoryTag;
  TAG       *emailTag;
  USER      *user;
  char      *string;

	string = malloc( strlen( baseName ) + 7 ); /* + ".rules" + '\0' */
  if( string == NULL )
    {
    fprintf( stderr, "Memory allocation failure in readUserList()\n" );
    return;
    }

  sprintf( string, "%s.list", baseName );
  file = newGroupFile( string );
  if( file == NULL )
    {
    printErrorMessage( "readUserList: Error",
      "newGroupFile() returned NULL in passwordCheck()" );
    /* fprintf( stderr, "newGroupFile() returned NULL in passwordCheck()\n" ); */
    free( string );
    return;
    }
  
  rules = newRules();
  if( rules == NULL )
    {
    printErrorMessage( "readUserList: Error", 
			"newRules() returned NULL in passwordCheck()" );
    /* fprintf( stderr, "newRules() returned NULL in passwordCheck()\n" ); */
    free( string );
    disposeGroupFile( file );
    return;
    }

  sprintf( string, "%s.rules", baseName );
  parseRuleFile( rules, string );
  /* printRules( rules ); */
  
  /* -- No error checking performed */

  group = getGroup( file, rules );
  while( group != NULL )
    {
    userTag = lookupTag( group, "user" );
    nameTag = lookupTag( group, "full-name" );
    passwordTag = lookupTag( group, "password" );
    homeDirectoryTag = lookupTag( group, "home-directory" );
    emailTag = lookupTag( group, "email" );

    if( userTag != NULL && nameTag != NULL && passwordTag != NULL && 
			  homeDirectoryTag != NULL && emailTag != NULL )
      {
      user = newUser( userTag->body, nameTag->body, passwordTag->body, 
        homeDirectoryTag->body, emailTag->body );
      addUserToUserList( user, userList );
      }

    discardGroup( file );
    group = getGroup( file, rules );
    }

  disposeGroupFile( file );
  disposeRules( rules );
  free( string );
  }


void printUserList( FILE *passwdFile, USERLIST *userList )
  {
	USERITEM *userItem;
	FILE *file;

	file = passwdFile;

	if (passwdFile == NULL)
		{
		file = stdout;
		}

	userItem = userList -> firstUserItem;
  while( userItem != NULL )
    {
		fprintf( file, "user: %s\n",      userItem -> user -> name );
		fprintf( file, "full-name: %s\n", userItem -> user -> fullName );
		fprintf( file, "email: %s\n",     userItem -> user -> email );
		fprintf( file, "password: %s\n",  userItem -> user -> password );
		fprintf( file, "home-directory: %s\n", userItem -> user -> homeDirectory );
		userItem = userItem -> nextUserItem;
		if (userItem != NULL) fprintf( file, "\n");
    }

  }


int addNewUser( char *name, char *fullName, char *password, 
							  char *homeDirectory, char *email, USERLIST *userList )
  {
	USER *user;

	user = newUser( name, fullName, password, homeDirectory, email );
	if (user == NULL)
		{
		fprintf( stderr, "addNewUser: error adding user\n");
		return -1;
		}

	addUserToUserList( user, userList );

	return 0;
	}


int setUserPassword( char *name, char *password, USERLIST *userList )
  {
  USER *user;

	if (userList == NULL)
		{
		fprintf(stderr,"setUserPassword: don''t have a userList\n");
		return -1;
		}

  user = lookupUser( name, userList );
	if (user == NULL)
		{
		fprintf(stderr,"setUserPassword: cannot find user \"%s\"\n", name);
		return -1;
		}

  user->password = strdup(password);
  stripTrailingSpace( user->password );

	return 0;
  }


void forwardUserTag( char *scriptname, USER *user, char *tag )
  {
    int len;
    char *format = "QUERY_STRING=user=%s&password=%s&%s";
    char *buff;

    len = strlen(format) + strlen(user->name) + strlen(user->password) + strlen(tag);
    buff = malloc(len);

    if (buff == NULL)
      {
	fprintf( stderr,"forwardUserTag: malloc failure\n" );
	return;
      }

    putenv( "REQUEST_METHOD=GET" );
    sprintf( buff, format, user->name, user->password, tag );
    putenv( buff );
    execl( scriptname, scriptname, NULL );
  }


void forwardUserDir( char *scriptname, USER *user, char *directory )
  {
    int len;
    char *format = "QUERY_STRING=user=%s&password=%s&directory=%s";
    char *buff;

    len = strlen(format) + strlen(user->name) + strlen(user->password) + strlen(directory);
    buff = malloc(len);

    if (buff == NULL)
      {
	fprintf( stderr,"forwardUserDir: malloc failure\n" );
	return;
      }

    putenv( "REQUEST_METHOD=GET" );
    sprintf( buff, format, user->name, user->password, directory );
    putenv( buff );
    execl( scriptname, scriptname, NULL );
  }


void forwardUser( char *scriptname, USER *user )
  {
#ifndef OLD_CODE
    int len;
    char *format = "QUERY_STRING=user=%s&password=%s";
    char *buff;

    len = strlen(format) + strlen(user->name) + strlen(user->password);
    buff = malloc(len);

    if (buff == NULL)
      {
	fprintf( stderr,"forwardUser: malloc failure\n" );
	return;
      }

    putenv( "REQUEST_METHOD=GET" );
    sprintf( buff, format, user->name, user->password );
    putenv( buff );
    execl( scriptname, scriptname, NULL );
#else
    fprintf( output, "Location: "
	     LIST_PROBLEM "?user=%s&password=%s\n\n", 
	     user->name, user->password );
#endif
  }

/* -- Public Methods ----------------------------------------------------------
*/

#if 0

void printJobHeading( FILE *output, USER *user )
  {
  char       *fileName;
  char       *jobDescription;
  GROUPFILE  *file;
  RULES      *rules;
  GROUP      *group;
  TAG        *tag;
  /* Constants */
  /*Now defined in Makefile*/
  /*const char *RuleFileName = "/usr/people/poor/CMISS/FormInterface/job.rules";*/
  /*const char *JobFileName  = "job.info";*/
  
  fileName = malloc( strlen( user->homeDirectory ) + strlen( JobFileName ) +
    2 );   /*  "/user/bob" + "/" + "job.info" + "\0" */
  if( fileName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in printJobHeading()\n" );
    return;
    }

  sprintf( fileName, "%s/%s", user->homeDirectory, JOB_FILE_NAME );

  file = newGroupFile( fileName );
  if( file == NULL )
    {  
    fprintf( stderr, "newGroupFile() reutned NULL in printJobHeading()\n" );
    free( fileName );
    return;
    }

  rules = newRules();
  if( rules == NULL )
    {
    fprintf( stderr, "newRules() returned NULL in printJobHeading()\n" );
    free( fileName );
    disposeGroupFile( file );
    return;
    }
  
  parseRuleFile( rules, RULE_FILE_NAME );
    /* -- NB. No error checking is preformed at this point */

  printf("c<p>");fflush(stdout);
  group = getGroup( file, rules );
  printf("c<p>");fflush(stdout);
  if( group == NULL )
    {
    fprintf( stderr, "getGroup() returned NULL in printJobHeading()\n" );
    free( fileName );
    disposeGroupFile( file );
    return;
    }

  tag = lookupTag( group, "job" );
  if( tag == NULL )
    {
    fprintf( stderr, "getTag() returned NULL in printJobHeading()\n" );
    free( fileName );
    disposeGroupFile( file );
    return;
    }

  fprintf( output, "<CENTER>" );
  fprintf( output, "<TABLE BGCOLOR=\"#ddbb88\" BORDER=\"0\" "
    "CELLPADDING=\"3\" CELLSPACING=\"0\"><TR><TD>" );
  fprintf( output, "<TABLE BGCOLOR=\"#ffffcc\" BORDER=\"0\" "
    "CELLPADDING=\"3\" CELLSPACING=\"0\">" );
  fprintf( output, "<TR><TD ALIGN=\"RIGHT\"><FONT SIZE=\"+1\" FACE=\"Helvetica\"><B>Username:</B></FONT></TD><TD><FONT SIZE=\"+1\" FACE=\"Helvetica\">%s</FONT></TD></TR>",
    user->name );
  fprintf( output, "<TR><TD ALIGN=\"RIGHT\"><FONT SIZE=\"+1\" FACE=\"Helvetica\"><B>Current Job:</FONT></B></TD><TD><FONT SIZE=\"+1\" FACE=\"Helvetica\">%s</FONT></TD></TR>",
    tag->body );
  fprintf( output, "</TABLE>" );
  fprintf( output, "</TD></TR></TABLE>" );
  fprintf( output, "</CENTER>" );

	fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">",user->name);
	fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">",user->password);


  free( fileName );
  disposeGroupFile( file );
  }

#endif 

