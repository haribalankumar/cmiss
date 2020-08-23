/*
 *  Password.c
 *
 *  Change passwds in the CMISS Internet Submission System
 */

/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <unistd.h>
#include <pwd.h>
#include "user-managment.h"
#include "transfer.h"
#include "webcmiss.h"
#include "password.h"

/* -- Program Constants -------------------------------------------------------
*/
#define SETUID

/* -- Module Private Methods --------------------------------------------------
*/

USERLIST *getPasswordList( FILE **passwordLock )
{
	USERLIST *userList;
	struct flock fl;

	fl.l_type   = F_WRLCK;
	fl.l_whence = SEEK_SET;
	fl.l_start  = 0;
	fl.l_len    = 0;

	/* Lock the password file */
  *passwordLock = fopen( ESU_MASTER_USER_LIST ".lock", "w+" );
	if (*passwordLock == NULL) {
    fprintf(stderr,"getPasswordList: cannot open password lock file \"%s\"\n", 
			ESU_MASTER_USER_LIST ".lock");
    return NULL;
  }
	if (fcntl(fileno(*passwordLock),F_SETLKW,&fl) < 0) {
    fprintf(stderr,"getPasswordList: cannot get lock on password file \"%s\" %d\n", 
			ESU_MASTER_USER_LIST ".lock", fileno(*passwordLock) );
		perror("fcntl");
    return NULL;
	}

  /* Open and read the password file */
  userList = newUserList();
  if( userList == NULL ) {
    fprintf( stderr, "getPasswordList: newUserList() returned NULL\n" );
    return NULL;
	}
	readUserList( userList, ESU_MASTER_USER_LIST );

  return userList;
}


int putPasswordList( USERLIST *userList, FILE *passwordLock )
{
	FILE *passwordNew;
	struct flock fl;

	fl.l_type   = F_UNLCK;
	fl.l_whence = SEEK_SET;
	fl.l_start  = 0;
	fl.l_len    = 0;

	passwordNew = fopen( ESU_MASTER_USER_LIST ".new", "w" );
  if (passwordNew == NULL) {
		fprintf(stderr, "putPasswordList: cannot open password file \"%s\"\n", 
			ESU_MASTER_USER_LIST ".new");
    return -1;
	}

	/* Write the password file */
	printUserList( passwordNew, userList );
	fclose(passwordNew);
	if (chmod(ESU_MASTER_USER_LIST ".new",S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)) {
		perror("putPasswordList: chmod");
		return -1;
	}

	/* Move the old password file to .back, and replace with the new file */
	if (move(ESU_MASTER_USER_LIST ".list",ESU_MASTER_USER_LIST ".back")) {
		perror( "putPasswordList: move[1]" );
    return -1;
	}
	if (move(ESU_MASTER_USER_LIST ".new",ESU_MASTER_USER_LIST ".list")) {
		perror( "putPasswordList: move[2]" );
    return -1;
	}

	/* Remove the password lock */
	fcntl(fileno(passwordLock),F_SETLKW,&fl);
	if (remove(ESU_MASTER_USER_LIST ".back"))	{
		perror( "putPasswordList: remove" );
    return -1;
	}

	return 0;
}



/* -- Public Methods ----------------------------------------------------------
*/

int changePassword( char *user, char *password )
{
	USERLIST *userList;
	FILE     *passwordLock;

	/* Open password file, get the user list */
	userList = getPasswordList( &passwordLock );
	if (userList == NULL) {
		fprintf(stderr,"changePassword: Error getting userList\n");
		return -1;
	}

	/* Make sure we have a username */
	if (user == NULL || password == NULL)	{
		fprintf(stderr,"changePassword: error in args\n");
		return -1;
	}

	/* Check username exists */
	if (lookupUser( user, userList ) == NULL) {
		fprintf(stderr,"changePassword: no such user\n");
		return -1;
	}

	/* Change the password */
	if (setUserPassword( user, password, userList )) {
		fprintf(stderr,"changePassword: error setting password\n");
		return -1;
	}

	/* Write out the new Password file */
	if (putPasswordList( userList, passwordLock )){
		fprintf(stderr,"changePassword: error writting password file\n");
		return -1;
	}

	return 0;
}


int addUser( char *user, char *password, char *email, char *fullName, char *home )
{
	USERLIST *userList;
	FILE     *passwordLock;

	/* Open password file, get the user list */
	userList = getPasswordList( &passwordLock );
	if (userList == NULL) {
		fprintf(stderr,"adduser: Error getting userList\n");
		return -1;
	}

	/* Make sure we have a username */
	if (user == NULL || password == NULL || email == NULL || fullName == NULL)	{
		fprintf(stderr,"adduser: error in args\n");
		return -1;
	}

	/* Check username exists */
	if (lookupUser( user, userList ) != NULL) {
		fprintf(stderr,"adduser: user \"%s\" already exists\n", user);
		return -1;
	}

	/* Add the user */
	if (addNewUser(user, fullName, password, home, email, userList)) {
		fprintf(stderr,"adduser: error ading user to list\n");
		return -1;
	}

	/* Write out the new Password file */
	if (putPasswordList( userList, passwordLock )){
		fprintf(stderr,"adduser: error writting password file\n");
		return -1;
	}

	return 0;
}

