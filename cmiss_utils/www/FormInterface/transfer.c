/*
 * Transfer.c
 *
 *    Routines to transfer files between the web server and hpc.
 */


/* -- Include Files -----------------------------------------------------------
*/

#undef SETUID

#include <sys/types.h>
#ifdef SETUID
#include <pwd.h>
    /* For: getpwnam() */
#endif

#include <stdio.h>
    /* For: fprintf(), fflush(), stdout, stderr, etc */

#include <stdlib.h>
    /* For: malloc(), free(), NULL */

#include <unistd.h>
    /* For: chdir() */

#include <sys/stat.h>
    /* For: stat() */

#include <sys/wait.h>
    /* For: wait() */

#include <fcntl.h>
    /* For: open() */

#include "user-managment.h"
    /* For: USER */

#include "job.h"
    /* For: JOB */

#include "host.h"
    /* For: HOST, HOSTITEM etc */

#include "webcmiss.h"
    /* For: SUBMISSION_SERVER, ABORT_SERVER etc */

#include "transfer.h"

#define BUFFER_SIZE 1024
#define TARFILE     ".transfer.tar"

#define UNDEBUG



int move(char *file1, char *file2)
{
	if (link(file1,file2))
		return -1;

	if (remove(file1))
		return -1;

	return 0;
}


int copy(char *file1, char *file2)
{
	FILE *fd1, *fd2;
	int c;

	fd1 = fopen(file1,"r");
	if (fd1 == NULL)
		return -1;
	fd2 = fopen(file2,"w");
	if (fd2 == NULL)
		return -1;

	while ((c = fgetc(fd1)) >= 0)
		fputc((char) c,fd2);

	fclose(fd1);
	fclose(fd2);

	return 0;
}


void checkDirectories(char *dir)
{
  struct stat file_info;

  if(stat(dir,&file_info) != 0)
  {
    if (mkdir(dir,S_IRWXU)) {
      fprintf(stderr,"checkDirectories: \"%s\" ", dir);
      perror("mkdir");
      exit(1);
    }
  }
}


static void receiveTarFile(FILE *sock,char *filename)
{
  FILE *fd;
  size_t readResult;
  char buf[BUFFER_SIZE];

  /* Create the tar file to write into */
  if((fd = fopen(filename, "wb")) == NULL)
  {
    fprintf(stderr,"receiveTarFile: ");
    perror("fopen");   
    exit(1);
  }

  /* Read from the standard input the tansferred tar file */
  while ((readResult = fread(buf, sizeof(char), BUFFER_SIZE, sock)) != 0)
  {
    if (fwrite(buf, sizeof(char), readResult, fd) != readResult)   
    {
      fprintf(stderr,"receiveTarFile: ");
      perror("fwrite");
      exit(1);
    }
  }
  fflush(fd);

  /* Close */
  if(fclose(fd)) 
  {
    fprintf(stderr,"receiveTarFile: ");
    perror("fclose"); 
    exit(1); 
  }
}


static void sendTarFile(FILE *sock, char *filename)
{
  FILE *fd;
  size_t readResult;
  char buf[BUFFER_SIZE];
  
  /* Open the file, and get a file handle */
  if((fd = fopen(filename,"rb")) == NULL)
  {
    fprintf(stderr,"sendTarFile: ");
    perror("fopen");
    exit(1);
  }

  /* Send the file */
  while ((readResult = fread(buf, sizeof(char), BUFFER_SIZE, fd)) != 0)
  {
    if (fwrite(buf, sizeof(char), readResult, sock) != readResult) 
    {
      fprintf(stderr,"sendTarFile: ");
      perror("fwrite");
      exit(1);
    }
  }
  fflush(sock);

  /* Close the file */
  if (fclose(fd)){
    fprintf(stderr,"sendTarFile: ");
    perror("fclose");
    exit(1);
  }
}


void receiveFiles(FILE *fd_in)
{
  /* Make sure we don't have a tar file to start with .... */
  system("/bin/rm -f " TARFILE);

  /* Receive the files */
  receiveTarFile(fd_in, TARFILE);

  /* Extract the input files */
  if (system("tar xf " TARFILE " 1>&2")) {
    fprintf(stderr,"receiveFiles: error untaring " TARFILE "\n");
    return;
  }

  /* Remove the tar file */
#ifndef DEBUG
  system("/bin/rm -f " TARFILE);
#endif
}


void sendFiles(FILE *fd_out)
{
  /* Make sure we don't have a tar file to start with .... */
  system("/bin/rm -f " TARFILE);

  /* Create the tar file */
  if (system("tar cf " TARFILE " *")) {
    fprintf(stderr,"sendFiles: error taring " TARFILE "\n");
    return;
  }

  /* Send the files */
  sendTarFile(fd_out, TARFILE);
  fflush(fd_out);

#ifndef DEBUG
  /* Remove the tar file */
  system("/bin/rm -f " TARFILE);
#endif
}


int submitJob(FILE *err, char *version, char *threads, 
	      HOST *host, JOB *job, USER *user)
{
  int pipefds_1[2], pipefds_2[2], pipe_read, pipe_write, ierr, len;
  int fd0, fd1, fd2;
  pid_t pid1, pid2;
  time_t startTime;
  char command[BUFFER_SIZE],line[BUFFER_SIZE],remote_pid[BUFFER_SIZE];
  char timestr[BUFFER_SIZE];
  FILE *fd_read, *fd_write;
  char *server;
  time_t start_time, stop_time;
  double wall_time, cpu_time;
  char *ctime_string;
  FILE *logfile;  
#ifdef SETUID
  struct passwd *pass;
#endif


  signal(SIGHUP,SIG_IGN);

  /* Construct our remote server name */
  len = strlen(host->homeDirectory) + 6 + strlen(SUBMIT_SERVER);
  server = malloc(len);
  if (server == NULL)
  {
    fprintf(stderr, "submitJob: malloc()\n");
    exit(0);
  }
  sprintf(server, "%s/bin/%s", host->homeDirectory, SUBMIT_SERVER );

  /* Open a pipe */
  if (pipe(pipefds_1) || pipe(pipefds_2)) {
    fprintf(err,"submitJob: pipe error\n");
    perror("pipe");
    return -1;
  }

  /* Run the job */
  if ((pid1 = fork()) == -1) {
    fprintf(err,"submitJob: forking error\n");
    perror("fork");
    return -1;
  }

  /* Child process */
  if (pid1 == 0) {

    /* Set stdin/stdout to the pipe file descriptors */
    close(0);
    dup(pipefds_1[0]);
    close(pipefds_1[0]);
    close(pipefds_1[1]);

    close(1);
    dup(pipefds_2[1]);
    close(pipefds_2[0]);
    close(pipefds_2[1]);

    /* close(2);*/

#ifdef SETUID
    /* SUID to webcmiss */
    pass = getpwnam("norris");
    setgid(pass->pw_gid);
    setuid(pass->pw_uid);
#endif

    /* Run the command */
    sprintf(command, "%s@%s", host->username, host->hostname);
    /* fprintf( stderr, "running \"%s %s %s\"\n", "rsh", command, server ); */
		
    /* if (ierr = execl("/usr/local/bin/ssh","ssh", command, server, NULL )) { */
    if (ierr = execl("/usr/bsd/rsh", "rsh", command, server, NULL )) {
      fprintf(stderr,"submitJob: command \"%s %s %s\" returns %d\n",
	      "rsh", command, server, ierr);
    }
    exit(0);
  }

  /* Parent process */
  close(pipefds_1[0]);
  close(pipefds_2[1]);
  
  pipe_read  = pipefds_2[0];
  pipe_write = pipefds_1[1];

  fd_read  = fdopen(pipe_read,  "r");
  fd_write = fdopen(pipe_write, "w");
  if (fd_read == NULL || fd_write == NULL)
  {
    fprintf(err, "submitJob: Can't create file pointer\n");
    return -1;
  }

  /* Get ACK/NAK response from the server */
  if (fgets(line, sizeof(line), fd_read) == NULL) 
  {  
    if (!feof(fd_read)) 
    { 
      fprintf(err, "submitJob: Error reading response\n"); 
      return -1;
    } 
  }

  switch (line[0])
  {
    case '+':
    {
      strcpy(remote_pid,&line[1]);
      break;
    }
    case '-':
    {
      fprintf(err,"Got NAK: %s\n", &line[1]);
      return -1;
    }
    default:
    {
      fprintf(err,"Got unknown response: \"%s\"\n", line);
      return  -1;
    }
  }

  /* Send user's name */
  if (fprintf(fd_write, "%s,%s,%s,%s,%s\n", user->name, user->email,
	      host->homeDirectory, version, threads) < 0)
  {
    perror("write");
    return -1;
  }
  fflush(fd_write);

  /* Check if a job is running */
  if (fgets(line, BUFFER_SIZE, fd_read) == NULL) 
  {
    perror("read[1]");
    return -1;
  }
  
  if (strcmp(line,"SENTINAL\n") == 0) {
    /* printf("submission: SENTINAL present -- no job will be run\n"); */
    return 1;
  }
  else if (strcmp(line,"ABSENT  \n") == 0) {
    /* printf("submission: no SENTINAL present -- job will be run\n"); */
  }
  else {
    /* printf("submission: cannot understand \"%s\"\n", line); */
    fprintf(err,"submission: cannot understand \"%s\"\n", line);
    return -1;
  }

  /* Send the input files */
  sendFiles(fd_write);

  /* Close the pipe */
  fclose(fd_write);

  /* Update the Job info file */
  startTime = time(NULL);
  sprintf(timestr,"%d",startTime);
  job->status = "running";
  job->host = host->hostname;
  job->pid = remote_pid;
  job->startTime = timestr;
  writeJobFile( user, job );

  /* Now get the returned files */
  fflush(stdin);
  fflush(stdout);
  if ((pid2 = fork()) == -1) {
    fprintf(err,"submitJob: forking error\n");
    perror("fork");
    return -1;
  }

  /* Child process */
  if (pid2 == 0) {

    /* Copy IO streams -- otherwise process either dies, or 
       tries to write to the html output */
    fd0 = open("/dev/null",O_RDONLY);
    close(0);
    dup(fd0);
    fd1 = open("/dev/null",O_WRONLY);
    close(1);
    dup(fd1);
    fd2 = open("/dev/null",O_WRONLY);
    close(2);
    dup(fd2);

    /* Get the run times */
    chdir("..");
    if (fgets(line, BUFFER_SIZE, fd_read) == NULL)
    {
      perror("read[2]");
      return -1;
    }
    sscanf(line, "%d,%d,%lf,%lf\n", &start_time, &stop_time, &wall_time, &cpu_time );
    ctime_string = ctime(&start_time);
    ctime_string[24] = '\0';

    /* Write them to the log file */
    sprintf( line, "%s/job.log", user->homeDirectory );
    logfile = fopen(line,"a");
    fprintf( logfile, "%s %12.3f %12.3f %8.8s %10.10s %4.4scpu %s\n", ctime_string,
      wall_time, cpu_time, host->hostname, version, threads, job->description );
    fclose( logfile );

    /* Return the files */
    chdir("Output/");
    system("rm -f * 2> /dev/null 1> /dev/null");
    receiveFiles(fd_read);
    sprintf(command,"cp `ls ../Input/%s.com|grep -v -e example_ "
	    "-e test_output.com` . 2> /dev/null 1> /dev/null", "*");
    system(command);
    system("chmod 644 * 2> /dev/null 1> /dev/null");

    /* Update the job file */
    startTime = time(NULL);
    sprintf(timestr,"%d",startTime);
    job->status = "finished";
    job->host = host->hostname;
    job->pid = "0";
    job->startTime = timestr;
    writeJobFile( user, job );

    exit(0);
  }

  return 0;
}

int abortJob(FILE *err, HOST *host, JOB *job, USER *user, int level)
{
  int pipefds_1[2], pipefds_2[2], pipe_read, pipe_write, ierr, len;
  pid_t pid;
  time_t startTime;
  char command[BUFFER_SIZE],line[BUFFER_SIZE],remote_pid[BUFFER_SIZE];
  char timestr[BUFFER_SIZE];
  FILE *fd_read, *fd_write;
  char *server;


  signal(SIGHUP,SIG_IGN);

  /* Construct our remote server name */
  len = strlen(host->homeDirectory) + 6 + strlen(ABORT_SERVER);
  server = malloc(len);
  if (server == NULL)
    {
      fprintf(stderr, "submitJob: malloc()\n");
      exit(0);
    }
  sprintf(server, "%s/bin/%s", host->homeDirectory, ABORT_SERVER );

  /* Open a pipe */
  if (pipe(pipefds_1) || pipe(pipefds_2)) {
    fprintf(err,"abortJob: pipe error\n");
    perror("pipe");
    return -1;
  }

  /* Run the job */
  if ((pid = fork()) == -1) {
    fprintf(err,"abortJob: forking error\n");
    perror("fork");
    return -1;
  }

  /* Child process */
  if (pid == 0) {

    /* Set stdin/stdout to the pipe file descriptors */
    close(0);
    dup(pipefds_1[0]);
    close(pipefds_1[0]);
    close(pipefds_1[1]);

    close(1);
    dup(pipefds_2[1]);
    close(pipefds_2[0]);
    close(pipefds_2[1]);

    /* close(2); */

    /* Run the command */
    sprintf(command, "%s@%s", host->username, host->hostname);
    /* fprintf( stderr, "running \"%s %s %s\"\n", "rsh", command, server ); */

    if (ierr = execl("/usr/bsd/rsh", "rsh", command, server, NULL )) {
      fprintf(stderr,"abortJob: command \"%s %s %s\" returns %d\n",
	      "rsh", command, server, ierr);
    }
    exit(0);
  }

  /* Parent process */
  close(pipefds_1[0]);
  close(pipefds_2[1]);
  
  pipe_read  = pipefds_2[0];
  pipe_write = pipefds_1[1];

  fd_read  = fdopen(pipe_read,  "r");
  fd_write = fdopen(pipe_write, "w");
  if (fd_read == NULL || fd_write == NULL)
  {
    fprintf(err, "abortJob: Can't create file pointer\n");
    return -1;
  }

  /* Get ACK/NAK response from the server */
  if (fgets(line, sizeof(line), fd_read) == NULL) 
  {  
    if (!feof(fd_read)) 
    { 
      fprintf(err, "abortJob: Error reading response\n"); 
      return -1;
    } 
  }

  switch (line[0])
  {
    case '+':
    {
      strcpy(remote_pid,&line[1]);
      break;
    }
    case '-':
    {
      fprintf(err,"Got NAK: %s\n", &line[1]);
      return -1;
    }
    default:
    {
      fprintf(err,"Got unknown response: %s\n", line);
      return  -1;
    }
  }

  /* Send user's name */
  if (fprintf(fd_write, "%s,%s,%s,%d\n", user->name, user->email,
	      host->homeDirectory, level) < 0)
  {
    perror("write");
    return -1;
  }
  fflush(fd_write);

  /* Update the job info */
  startTime = time(NULL);
  sprintf(timestr,"%d",startTime);
  job->status = "finished";
  job->host = host->hostname;
  job->pid = "0";
  job->startTime = timestr;
  writeJobFile( user, job );

  return 0;
}
