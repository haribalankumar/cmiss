/*******************************************************************************
FILE : abort-server.c

LAST MODIFIED : 30 April 2002

DESCRIPTION :
This program reads the standard input the users name and email address.
It then tries to abort the previous CMISS job started by the user.
==============================================================================*/
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <errno.h>
#include <unistd.h>
#include <netinet/in.h>
#include <limits.h> 
#include <netdb.h> 
#include <arpa/inet.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <pwd.h>
     /* For: getpwnam() */

#include <sys/utsname.h>
     /* For: uname() */

#include <grp.h>
#include <fcntl.h>
#include "transfer.h"
#include "webcmiss.h"

#define BUFFER_SIZE 1024
#undef SETUID

extern int errno;
     /* For: error number */


static void shiftError(char *err_file, char *outdir)
{
  char command[BUFFER_SIZE];
  
  sprintf(command,"/bin/mv %s %s/runtime_stderr_output", err_file, outdir); 
  system(command);
}


int main(int argc, char **argv)
{
  char email[BUFFER_SIZE], path[BUFFER_SIZE], user[BUFFER_SIZE];
  char indir[BUFFER_SIZE], outdir[BUFFER_SIZE], topdir[BUFFER_SIZE];
  char err_file[BUFFER_SIZE], buf[BUFFER_SIZE], buff[BUFFER_SIZE];
  FILE *sentinal;
  struct stat file_info;
  pid_t process_id, parent_id, child_id;
#ifdef SETDUID
  struct passwd *pass;
#endif
	mode_t mask_mode = S_IRWXG | S_IRWXO;
	struct utsname uname_data;
  int level;


	uname(&uname_data);

	/* Set umask: we want to keep to ourselves */
	(void) umask(mask_mode);

  /* Send acknowledgement */
  process_id = getpid();
  printf("+%d\r\n", process_id);
  fflush(stdout);

  /* Read the user's name and email address */
  if (fgets(buf, BUFFER_SIZE, stdin) == NULL) {  
    fprintf(stderr, "Error reading response\n"); 
    return -1;
  }
  sscanf(buf, "%[^,],%[^,],%[^,],%[^\n]", user, email, path, buff);
  level = atoi(buff);

  /* Open a file for the stderr stream */
  sprintf(err_file, "%s/people/stderr_abort_%d", path, process_id);
  if (freopen(err_file, "wa", stderr) != stderr) {
    perror("freopen");
  }

  fprintf(stderr,"\nThis is the standard error output of the of\n");
  fprintf(stderr,"job cancelling program on %s\n\n", uname_data.nodename);
  if (level) {
    fprintf(stderr,"Running at level %d (kill procs and clear files)\n\n", level);
  }
  else {
    fprintf(stderr,"Running at level %d (abort job)\n\n", level);
  }
  fprintf(stderr,"-------------------------------------------\n\n");
  fflush(stderr);
  
#ifdef SETUID
  pass = getpwnam(user);
  setgid(pass->pw_gid);
  setuid(pass->pw_uid);
#endif

  /* Check if the user has a directory, and input and output directories */
  sprintf(topdir, "%s/people/%s", path, user); 
  checkDirectories(topdir);
  sprintf(indir,  "%s/people/%s/Input", path, user); 
  checkDirectories(indir);
  sprintf(outdir, "%s/people/%s/Output", path, user); 
  checkDirectories(outdir);

  shiftError(err_file, outdir);

  /* Make the Input directory the current working directory */
  if(chdir(indir)) {
    perror("chdir");
    return -1;
  }

  /* Check if a SENTINAL file exists. If so we are running a job which
     we want to kill. Otherwise we just delete any files that might be
		 hanging around. */
  sprintf(topdir, "%s/SENTINAL", indir); 
  if(stat(topdir, &file_info) == 0)	{
		sentinal = fopen("SENTINAL", "r");
    fscanf(sentinal, "%d", &parent_id );
    fscanf(sentinal, "%d", &child_id );
		fclose(sentinal);
		sleep(1);
  }

  /* Kill the processes: we first kill the child. If we are running at level 1
     or above we also kill the parent, and delete all files. */
  kill(child_id, SIGKILL);
  if (!level)
  {
    return 0;
  }

  /* We are on a search and destry mission -- get rid of the parent and all files */
  sleep(2);
  kill(parent_id, SIGKILL);

  /* Delete existing files, and receive the new ones */
	if (chdir(indir)) {
		fprintf(stderr, "abort-server can't cd to \"%s\"\n", indir);
		return -1;
  }
	system("/bin/rm -f *");

	if (chdir(outdir)) {
		fprintf(stderr, "abort-server: can't cd to \"%s\"\n", outdir);
		return -1;
  }
	system("/bin/rm -f *");

  return 0;
}
