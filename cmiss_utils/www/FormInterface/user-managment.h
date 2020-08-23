/*
 *  User Managment.h
 *
 *    Public Interface to the User Managment routines
 */
#ifndef USER_MANAGMENT_H
#define USER_MANAGMENT_H


/* -- Required Include Files --------------------------------------------------
*/

#include <time.h>
    /* For: time() */

#include "groups.h"
    /* For: GROUP, GROUPFILE, setGroupFilePos(), etc, etc */

#include "rules.h"
    /* For: RULES, newRules(), etc */

#include "rule-parse.h"
    /* For: parseRules() */



/* -- Module Datatypes --------------------------------------------------------
*/

typedef struct
  {
  char *name;
  char *fullName;
  char *password;
  char *homeDirectory;
  char *email;
  }
USER;


typedef struct USERITEM
  {
  USER *user;

  struct USERITEM *nextUserItem;
  }
USERITEM;


typedef struct 
  {
  USERITEM *firstUserItem;
  }
USERLIST;



/* -- Module Constructor and Destructors --------------------------------------
*/

USER *newUser( char *name, char *fullName, char *password, char *homeDirectory, 
  char *email );
USER *duplicateUser( USER *user );
void disposeUser( USER *user );

USERITEM *newUserItem( USER *user );
void disposeUserItem( USERITEM *userItem );
void disposeUserItems( USERITEM *userItem );

USERLIST *newUserList( void );
void disposeUserList( USERLIST *userList );



/* -- Public Methods ----------------------------------------------------------
*/

USER *lookupUser( char *name, USERLIST *userList );
void addUserToUserList( USER *user, USERLIST *userList );
void readUserList( USERLIST *userList, char *baseName );
void printUserList( FILE *passwdFile, USERLIST *userList );
int addNewUser( char *name, char *fullName, char *password, 
							  char *homeDirectory, char *email, USERLIST *userList );
int setUserPassword( char *name, char *password, USERLIST *userList );
void forwardUser( char *scriptname, USER *user );
void forwardUserDir( char *scriptname, USER *user, char *directory );
void forwardUserTag( char *scriptname, USER *user, char *tag );

/*
void printJobHeading( FILE *output, USER *user );
*/


#endif
