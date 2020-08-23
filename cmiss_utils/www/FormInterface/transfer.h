#ifndef TRANSFER_H
#define TRANSFER_H


/* -- Required Include Files --------------------------------------------------
*/

#include "user-managment.h"
    /* For: USER */

#include "job.h"
    /* For: JOB */

#include "host.h"
    /* For: HOST, HOSTITEM etc */

void checkDirectories(char *dir);
void receiveFiles(FILE *fd_in);
void sendFiles(FILE *fd_out);

int move(char *file1, char *file2);
int copy(char *file1, char *file2);

int returnFiles(char *host, char *ruser, char *user, char *server);
int submitJob(FILE *err, char *version, char *threads, HOST *host, JOB *job, USER *user);
int abortJob(FILE *err, HOST *host, JOB *job, USER *user, int level);

#endif
