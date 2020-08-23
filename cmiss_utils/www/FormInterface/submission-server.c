/*******************************************************************************
FILE : submission-server.c

LAST MODIFIED : 01 May 2002

DESCRIPTION :
This program reads the standard input the users name and email address,
and the tared input files. The input files untared into the input directory.
After running CMISS the output files are copied to the output
directory, tarred and sent back to the webserver.

It is envisaged that in future the user will be able run multiple jobs,
This will require the input and output will go to an appropriate
directories (perhaps identified by process number). Work on this has begun
but is now commented out because it involves significant changes to
other code.
==============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/time.h>


#include <sys/wait.h>
     /* For: wait() */

#include <sys/utsname.h>
     /* For: uname() */

#include <unistd.h>
     /* For: fork(), pipe() etc */

#include <signal.h>
#include <errno.h>
#include <unistd.h>
#include <netinet/in.h>
#include <limits.h> 
#include <netdb.h> 
#include <arpa/inet.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pwd.h>
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
  
  sprintf(command, "/bin/mv %s %s/runtime_stderr_output", err_file, outdir); 
  system(command);
}


static void exitClean(int ierr, char *err_file, char *outdir)
{
  shiftError(err_file, outdir);

  exit(ierr);
}


int main(int argc, char **argv)
{
  char command[BUFFER_SIZE], email[BUFFER_SIZE], topdir[BUFFER_SIZE];
  char threads[BUFFER_SIZE], user[BUFFER_SIZE], path[BUFFER_SIZE];
  char version[BUFFER_SIZE], err_file[BUFFER_SIZE], indir[BUFFER_SIZE];
  char outdir[BUFFER_SIZE], buf[BUFFER_SIZE];
  FILE         *fd;
  FILE         *fp;
  FILE         *sentinal;
  int          statid;
  int          pipefds[2];
  char         *comfile = "cm_test.com";
  char         *time_info;
  struct stat file_info;
  pid_t process_id, cm_pid;
  time_t time_start, time_stop;
  double cpu_time, parent_time, child_time, wall_time;
  struct tms cpu_times;
#ifdef SETDUID
  struct passwd *pass;
#endif
  mode_t mask_mode = S_IRWXG | S_IRWXO;
  struct utsname uname_data;


  uname(&uname_data);

  /* Set umask: we want to keep to ourselves */
  (void) umask(mask_mode);

  /* Send acknowledgement */
  process_id = getpid();
  printf("+%d\r\n", process_id);
  fflush(stdout);

  /* Read the user's name and email address */
  if (fgets(buf, BUFFER_SIZE, stdin) == NULL)
  {
    fprintf(stdout, "submission-server: Error reading response\n"); 
    return 1;
  }
  sscanf(buf, "%[^,],%[^,],%[^,],%[^,],%[^\n]", user, email, path,
         version, threads);

  /* Open a file for the stderr stream */
  sprintf(err_file,"%s/people/stderr_submitter_%d", path, process_id);
  if (freopen(err_file, "w", stderr) != stderr) {
    perror("freopen");
  }
  
  fprintf(stderr, "\nThis is the standard error output of the of\n");
  fprintf(stderr, "job submission program on %s running %s\n\n",
          uname_data.nodename, uname_data.sysname);
  fprintf(stderr, "-------------------------------------------\n\n");
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

  /* Make the Input directory the current working directory */
  if(chdir(indir))
  {
    perror("chdir");
    exitClean(1, err_file, outdir);
  }

  /* Check if a SENTINAL file exists. If so we are running a job.
     Then we exit. Otherwise create a sentinel. */
  sprintf(topdir, "%s/SENTINAL", indir);
  if(stat(topdir, &file_info) == 0)
  {
    printf("SENTINAL\n");
    fflush(stdout);
    exitClean(0, err_file, outdir);
  }
  else {
    errno = 0;
    sentinal = fopen("SENTINAL", "w");
    fprintf(sentinal, "%d\n", process_id);
    fflush(sentinal);

    printf("ABSENT  \n");
    fflush(stdout);
  }

  /* If no sentinel, receive some tar files */

  /* Delete existing files, and receive the new ones */
  system("rm -f `ls | grep -v SENTINAL`");
  sprintf(command, "rm -f %s/%s", outdir, "*");
  system(command);
  receiveFiles(stdin);
  system("chmod u=rw *");

  /* Shift down the error file */
  shiftError(err_file, outdir);

  /* Create the new com file */
  sprintf(command,"( echo \"set fatal off;\" ; echo \"set num_threads %s;\" ;"
    " cat example_*.com test_output.com ) > %s", threads, comfile);   
  system(command);

  /* Touch all the files. We want them to be dated at the same time */
  /* as the file we create called "CURRENTDATE". This file is used in */
  /* a "find -newer" operation later on to seperate files with */
  /* modification dates newer than the input file (i.e. any output */
  /* files generated by the code run). */
  system("touch *"); 
  system("touch CURRENTDATE");

  /*  Sleep for two seconds because the machine may output files faster than */
  /*  the minimum clock resolution of the file system (and thus we may get */
  /*  new files with the same creation date as the CURRENTDATE file). */
  sleep(2);
  time_start = time(NULL);

  /* Create a pipe so we can talk to our code (it might get lonely) */
  if (pipe(pipefds) < 0) {
    perror("pipe");
    return -1;
  }
  
  /*  Run CMISS with the command file, and pipe any screen IO to job.out */
  /*  The epath support is in case people forget to remove the */
  /*  ;DOC lines from their comfiles.  */
  if ((cm_pid = fork()) == -1) {
    perror("fork");
    return -1;
  }

  /* Child process */
  if (cm_pid == 0) {

    /*
     * Make the read side of the pipe our standard input. Close
     * the write side of the pipe.
     */
    close(0);
    dup(pipefds[0]);
    close(pipefds[0]);
    close(pipefds[1]);

    /* Set the enviroment for the number of threads: no need to do
       this -- the number of threads is set in the com file above */
    sprintf(command, "NUM_THREADS=%s", threads);
    if (putenv(command)) {
      perror("putenv");
      return -1;
    }

    /* Set the stdout to job.out */
    if (freopen("job.out", "w", stdout) != stdout) {
      perror("freopen");
      return -1;
    }
    /* Pipe the stderr to job.out */
    close(2);
    dup(1);
    fclose(sentinal);

    /* exec */
    sprintf(command, "%s/%s/%s", path, "bin", version);
#ifdef DEBUG
    fprintf(stderr, "running %s/%s/%s -epath %s %s", path, "bin", version, indir, comfile);
#endif

    /* By doing this we ensure that the .profile file is read and ulimits are set correctly */
#ifdef OLD
    if (execl(command, version, "-epath", indir, comfile, NULL)) {
#else
    if (execl("/bin/ksh", "-ksh", "exec", "time", command, "-epath", indir, comfile, NULL)) {
#endif
      perror("execl");
      return -1;
    }
    return 0;
  }
  
  /* Write the cm process id to file */
  fprintf(sentinal, "%d\n", cm_pid);
  fclose(sentinal);

  /* Write a calming "quit" command to the running process */
  close(pipefds[0]);
  fp = fdopen(pipefds[1], "w");
  fprintf(fp, "\nquit\n");
  fclose(fp);

  /* Sit and twiddle our thumbs until the child finishes */
  wait(&statid);

  /*  Copy the output files to the Output directory */
  sprintf(command, "/bin/cp `find ./ -type f -newer CURRENTDATE -print"
          " | grep -v -e SENTINAL -e ./tmp. ` %s", outdir); 
  system(command);

  /* Send back the process times */
  times(&cpu_times);
  parent_time = ( (double) ( cpu_times.tms_utime + cpu_times.tms_stime ) )
              / (double) sysconf(_SC_CLK_TCK) ;
  child_time  = ( (double) ( cpu_times.tms_cutime + cpu_times.tms_cstime ) )
              / (double) sysconf(_SC_CLK_TCK) ;
  cpu_time  = parent_time + child_time;
  time_stop = time(NULL);
  wall_time = (double) ( time_stop - time_start );
  printf("%d,%d,%f,%f\n", time_start, time_stop, wall_time, cpu_time);
  
  /* Make the Output directory the current working directory */
  if(chdir(outdir))
  {
    perror("chdir");
    return 1;
  }

  /* Send the output files to RETURN_HOST */
  sendFiles(stdout);
  fclose(stdout);

  /* Mail the user notification of job completion */
  sprintf(command, "%s/bin/mail-user.sh %s %s", path, user, email );
  system(command);

  /* Delete Sentinal and files ready for another job */
  sprintf(command, "/bin/rm -f %s/%s", outdir, "*");
  system(command); 
  sprintf(command, "/bin/rm -f `ls %s/%s | grep -v SENTINAL`", indir, "*" );
  system(command);
  sprintf(command, "/bin/rm -f %s/SENTINAL", indir );
  system(command);

  /* Get the cpu times and the time of day */
  time_info = ctime(&time_stop);
  time_info[24] = '\0';

  /* Log the job */
  sprintf(topdir, "%s/submission-server.log", path);
  if((fd = fopen(topdir,"a")) == NULL) {
    perror("fopen");
    return 1;
  }

  fprintf(fd, "%s %s %s/%s\t%.3f %.1f\n", time_info, user, 
          version, threads, cpu_time, wall_time );
  if(fclose(fd)) {
    perror("fclose");
    return 1;
  }

  return 0;
}
