/*
 *  webcmiss_passwd.c
 *
 *  Change passwds/Add users to the CMISS Internet Submission System
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>
#include <crypt.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <grp.h>
#include <pwd.h>
#include "user-managment.h"
#include "transfer.h"
#include "password.h"
#include "webcmiss.h"


/* -- Program Constants -------------------------------------------------------
*/
#define LINESIZE 200


char *getNewPass(char *user)
{
	char *pass, *pass1, *pass2;

	/* Get password -- ask twice, and make sure they match */
	if ( ((pass = getpass(" Enter password: ")) == NULL) ||
			 ((pass1 = strdup(pass)) == NULL)) {
		perror("getpass");
		return NULL;
	}

	if (strlen(pass1) < 6) {
		fprintf(stderr," password should be at least 6 characters long.\n");
		return NULL;
	}

	if ( ((pass = getpass(" Repeat password: ")) == NULL) ||
			 ((pass2 = strdup(pass)) == NULL)) {
		perror("getpass");
		return NULL;
	}

	if (strcmp(pass1,pass2) != 0)	{
		fprintf(stderr," passwords do not match.\n");
		return NULL;
	}
	pass = crypt(pass1,user);

	return pass;
}


void createFiles(char *user)
{
	FILE *jobFile;
	char *file;
	mode_t mask_mode = S_IWGRP | S_IWOTH;
	mode_t dir_mode  = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
	struct passwd *nobody_user;
	struct group  *nobody_group;

	file = malloc( sizeof(ESU_WEB) + sizeof("%speople/%s/job.info")
							 + sizeof(user));
	if (file == NULL)	{
		perror("malloc");
		exit(-1);
	}

	/* Find nobody */
	if ((nobody_user = getpwnam("nobody")) == NULL)	{
    fprintf( stderr, "getpwnam(): can't find group\n" );
    exit(-1);
	}
	if ((nobody_group = getgrnam("nobody")) == NULL) {
    fprintf( stderr, "getgrnam(): can't find group\n" );
    exit(-1);
	}

	/* Make the account directories and the job.info file */
	(void) umask(mask_mode);

  sprintf(file,"%speople/%s",ESU_WEB,user);
  if (mkdir(file,dir_mode))	{
		perror("mkdir");
		exit(-1);
	}
  sprintf(file,"%speople/%s/Input",ESU_WEB,user);
  if (mkdir(file,dir_mode))	{
		perror("mkdir");
		exit(-1);
	}
  if (chown(file,nobody_user -> pw_uid,nobody_group -> gr_gid))	{
		perror("chown");
		exit(-1);
	}

  sprintf(file,"%speople/%s/Output",ESU_WEB,user);
  if (mkdir(file,dir_mode))	{
		perror("mkdir");
		exit(-1);
	}
  if (chown(file,nobody_user -> pw_uid,nobody_group -> gr_gid))	{
		perror("chown");
		exit(-1);
	}

  sprintf(file,"%speople/%s/job.info",ESU_WEB,user);
  if ((jobFile = fopen(file,"w")) == NULL) {
		perror("fopen");
		exit(-1);
	}
	fprintf(jobFile, "Job:\n");
	fprintf(jobFile, "Original-Files:\n");
	fprintf(jobFile, "Status: finished\n");
	fprintf(jobFile, "Host: -\n");
	fprintf(jobFile, "PID: 0\n");
	fprintf(jobFile, "StartTime: 0\n");
	fclose(jobFile);
	if (chmod(file,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)) {
		perror("chmod");
		exit(-1);
	}
  if (chown(file,nobody_user -> pw_uid,nobody_group -> gr_gid)) {
		perror("chown");
		exit(-1);
	}

	sprintf(file,"%speople/%s/job.log",ESU_WEB,user);
  if ((jobFile = fopen(file,"w")) == NULL) {
		perror("fopen");
		exit(-1);
	}
	fprintf(jobFile, "# Start time                Wall time     CPU time"
    "     Host    Version    Proc Title\n");
	fclose(jobFile);
	if (chmod(file,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)) {
		perror("chmod");
		exit(-1);
	}
  if (chown(file,nobody_user -> pw_uid,nobody_group -> gr_gid)) {
		perror("chown");
		exit(-1);
	}

	free(file);

	return;
}


void printUsage(int ierr)
{
	fprintf(stderr,"webcmiss_passwd -c          create user\n"
                 "                -u <user>   change password for \"user\"\n" );
	exit(ierr);
}


/* -- Program Entry Point -----------------------------------------------------
*/

void main( int argc, char *argv[] )
{
#ifdef CHECK_USERNAME
  USERLIST  *userList;
#endif
  char line[LINESIZE];
	char *pass;
	int  i;
	enum operation_type { CREATE, CHANGE } operation = CHANGE;
	char *user  = NULL;
	char *name  = NULL;
	char *email = NULL;


	/* Set signals, in case we are interupted */
	signal(SIGHUP,exit);
	signal(SIGINT,exit);
	signal(SIGQUIT,exit);
	signal(SIGTERM,exit);


	/* Parse the args -- see what we are actually up to */
	for (i = 1; i < argc; i++) {
		if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-h") == 0))	{
			printUsage(0);
			exit(0);
		}
		else if (strcmp(argv[i],"-c") == 0) {
			operation = CREATE;
		}
		else if (strcmp(argv[i],"-u") == 0)	{
			i++;
			if (i >= argc) printUsage(-1);
			user = argv[i];
		}
		else if (strcmp(argv[i],"-e") == 0)	{
			i++;
			if (i >= argc) printUsage(-1);
			email = argv[i];
		}
		else if (strcmp(argv[i],"-n") == 0) {
			i++;
			if (i >= argc) printUsage(-1);
			name = argv[i];
		}
		else {
			printUsage(-1);
			exit(-1);
		}
	}

#ifdef CHECK_USERNAME
  /* Open and read the password file */
  userList = newUserList();
  if( userList == NULL ) {
    fprintf( stderr, "getPasswordList: newUserList() returned NULL\n" );
    exit(-1);
	}
	readUserList( userList, ESU_MASTER_USER_LIST );
#endif

	/* Two types of operation -- change password, and add user */

	/* Change password */
	if (operation == CHANGE) {

		/* Make sure we have a username */
		if (user == NULL)	{
			printf(" Enter the login name of the user: \n");
			gets(line);
			user = strdup(line);
		}

		/* Check username exists */
#ifdef CHECK_USERNAME
		if (lookupUser( user, userList ) == NULL) {
			fprintf(stderr,"webcmiss_passwd: no such user \"%s\".\n", user);
			exit(-1);			
		}
#endif

		/* Get password -- ask twice, and make sure they match */
		pass = getNewPass(user);
		if (pass == NULL) {
			exit(-1);
		}

		/* Change the password */
		if (changePassword( user, pass )) {
			fprintf(stderr,"webcmiss_passwd: error changeing password.\n");
			exit(-1);
		}
	}

  /* Add user */
  else {

		/* Prompt for incomplete information */
		printf(" CMISS Internet Submission System User Manager\n\n");
		if (user == NULL)	{
			printf(" Enter the login name of the new user: \n");
			gets(line);
			user = strdup(line);
		}

		/* Check username not yet used */
#ifdef CHECK_USERNAME
		if (lookupUser( user, userList ) != NULL) {
			fprintf(stderr,"webcmiss_passwd: user \"%s\" already exists.\n", user);
			exit(-1);			
		}
#endif

		if (name == NULL)	{
			printf(" Enter the full name of the new user: \n");
			gets(line);
			name = strdup(line);
		}

		if (email == NULL) {
			printf(" Enter email address: \n");
			gets(line);
			email = strdup(line);
		}

		/* Get password -- ask twice, and make sure they match */
		pass = getNewPass(user);
		if (pass == NULL) {
			exit(-1);
		}

    /* Add the account details to the password file */
		sprintf(line, "%speople/%s\n",ESU_WEB,user);
		if (addUser( user, pass, email, name, line )) {
			fprintf(stderr,"webcmiss: error in addUser.\n");
			exit(-1);
		}

		/* Make the account directories and the job.info file */
		createFiles(user);
	}

	exit(0);
}
