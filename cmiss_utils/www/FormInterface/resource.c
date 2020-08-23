/*
 *  Resource.c
 *
 *    Module for calculating user resources, both disk space and
 *    CPU time used.
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <sys/types.h>
    /* For: *_t types */

#include <stdio.h>
    /* For: printf(), etc */

#include <stdlib.h>
    /* For: malloc() */

#include <sys/stat.h>
    /* For: stat(), lstat(), S_ISDIR() etc */

#include <unistd.h>
    /* For: chdir(), fchdir(), etc */

#include <sys/dir.h>
    /* For: opendir(), etc [POSIX style] */

#include <fcntl.h>
    /* For: O_RDONLY, etc [used with open()] */

#include "host.h"
    /* For: HOSTLIST, etc */

#include "resource.h"


#define BUFFER_SIZE 512


/* -- Local methods -----------------------------------------------------------
*/

static long long getFileUsage(char *fileName)
{
	long long usage;
	struct stat fileStat;

	/* stat() the file. We use lstat() so as not to follow symbolic links */
	if (lstat(fileName,&fileStat)) {
		perror(fileName);
		return -1LL;
	}

	/* Get the size of the file */
	/* usage = fileStat.st_size; */
	/* Get the number of blocks used */
	usage = fileStat.st_blocks;

	return usage;
}


static long long getDirectoryUsage(char *dirName)
{
	DIR *dir;
	struct direct *dirEnt;
	struct stat dirStat, fileStat;
	int cwd;
	long long usage = 0LL;


	/* stat() the file. We use lstat() so as not to follow symbolic links */
	if (lstat(dirName,&dirStat)) {
		perror(dirName);
		return -1LL;
	}

	/* Check we have a directory to check out */
	if (!S_ISDIR(dirStat.st_mode)) {
		return 0LL;
	}

	/* Save the location of the current directory */
	if ((cwd = open(".",O_RDONLY)) <= 0) {
		perror("opencwd");
		return -1LL;
	}

	/* Change to the named directory */
	if (chdir(dirName)) {
		perror("chdir()");
		return 0LL;
	}

	/* Open the directory */
	if ((dir = opendir(".")) == NULL) {
		perror("opendir");
		return 0LL;
	}

	/* Loop over the directory entries */
	while ((dirEnt = readdir(dir)) != NULL) {

		/* Ignore the parent directory */
		if (strcmp(dirEnt->d_name,"..") != 0) {

			/* Stat the entry. We use lstat() so as not to follow symbolic links */
			if (!lstat(dirEnt->d_name,&fileStat)) {

				/* If a directory, recurse */
				if (S_ISDIR(fileStat.st_mode) && strcmp(dirEnt->d_name,".") != 0) {
					usage += getDirectoryUsage(dirEnt->d_name);
				}
				/* Else, get the size */
				else {
					usage += fileStat.st_blocks;
				}
			}
		}
	}

	/* Close the directory and clean up */
	closedir(dir);

	/* Return to the original directory */
	if (fchdir(cwd)) {
		perror("fchdir");
	}
	close(cwd);

	return usage;
}


/* -- Public interface -----------------------------------------------------------------
*/

long long getDiskUsage(char *nodeName)
{
	long long usage = 0LL;
	struct stat nodeStat;

	/* Don't want top open ".." */
	if (strcmp(nodeName,"..") == 0) {
		return 0LL;
	}

	/* stat() the file. We use lstat() so as not to follow symbolic links */
	if (lstat(nodeName,&nodeStat)) {
		perror(nodeName);
		return -1LL;
	}

	/* Check we have a file or a directory */
	if (S_ISDIR(nodeStat.st_mode)) {
		usage = getDirectoryUsage(nodeName);
	}
	else /* if (S_ISFILE(nodeStat.st_mode)) */ {
		usage = getFileUsage(nodeName);
	}

	return usage;
}


void fprintDiskUsage(FILE *file, long long diskUsage)
{
#define K ( 2 )
#define M ( K * 1024 )
#define G ( M * 1024 )
  double du;

  if ((du = (double) diskUsage / (double) G) > 1.0) {
    fprintf( file, "%.2f GB", du );
  }
  else if ((du = (double) diskUsage / (double) M) > 1.0) {
    fprintf( file, "%.2f MB", du );
  }
  else {
    fprintf( file, "%d kB", (int) diskUsage / 2 );
  }
}


int getCPUUsageByHost(char *filename, char *hostname,
											double *wall_usage, double *cpu_usage)
{
	FILE *file;
	char buff[BUFFER_SIZE], host[BUFFER_SIZE];
	double a, b;

	*wall_usage = 0.0;
	*cpu_usage = 0.0;

	if ((file = fopen( filename, "r" )) == NULL) {
		perror("logfile");
		return -1;
	}

	while (fgets(buff, BUFFER_SIZE, file))
	{
		if (buff[0] != '#') {
			sscanf( &(buff[25]), "%lf %lf %s", &a, &b, host );
			if (strcmp(host,hostname) == 0) {
				*wall_usage += a;
				*cpu_usage  += b;
			}
		}
	}

	return 0;
}


USAGELIST *allocUsageList(int nhost)
{
	USAGELIST *usageList;

	if ((usageList = calloc(1,sizeof(USAGELIST))) == NULL) {
		perror("calloc");
		return NULL;
	}

	usageList->nhost = nhost;

	if ((usageList->usage = calloc(nhost,sizeof(USAGE))) == NULL) {
		perror("calloc");
		return NULL;
	}

	return usageList;
}


void freeUsageList(USAGELIST *usageList)
{
	free(usageList->usage);
	free(usageList);
}


USAGELIST *getUsageFromLog(char *filename, HOSTLIST *hostList)
{
	USAGELIST *usageList;
	HOSTITEM  *hostItem;
	FILE      *file;
	char      buff[BUFFER_SIZE], host[9];
	int       i, nhost;
	double    wall, cpu;

	/* Get the number of hosts */
	nhost = 0;
	hostItem = hostList->firstHostItem;
  while( hostItem != NULL ) {
		nhost ++;
		hostItem = hostItem->nextHostItem;
	}

	/* Allocate memory to the usage struct, set hostnames, clear other datafields */
	usageList = allocUsageList(nhost);
	hostItem = hostList->firstHostItem;
	for (i=0;i<nhost;i++) {
		strcpy(usageList->usage[i].hostname,hostItem->host->hostname);
		usageList->usage[i].wall_total = 0.0;
		usageList->usage[i].cpu_total  = 0.0;

		usageList->usage[i].wall_time = 0.0;
		usageList->usage[i].cpu_time  = 0.0;
		usageList->usage[i].start_time[0]  = '\0';
		usageList->usage[i].version[0]     = '\0';
		usageList->usage[i].threads[0]     = '\0';
		usageList->usage[i].description[0] = '\0';

		usageList->usage[i].start_time[24] = '\0';
		usageList->usage[i].version[10]    = '\0';
		usageList->usage[i].threads[8]     = '\0';

		hostItem = hostItem->nextHostItem;
	}

	/* Open the log file */
	if ((file = fopen( filename, "r" )) == NULL) {
		perror("logfile");
		return NULL;
	}

	/* Loop over logfile entries, ignoreing comments */
	while (fgets(buff, BUFFER_SIZE, file))
	{

		/* Not a comment -- find the hostname */
		if (buff[0] != '#') {

			/* Read the host name, and work out which host we are in */
			sscanf( &(buff[50]), "%s", host );
			for (i=0;i<nhost;i++) {
				if (strcmp(host,usageList->usage[i].hostname) == 0) {
					break;
				}
			}
			/* Check we found a host: if not, ignore the entry */
			if (i >= nhost) continue;

			/* Read the entries from the log file. Sum the cpu/wallclock times */
			strncpy(usageList->usage[i].start_time,buff,24);
			sscanf(&(buff[24]),"%lf %lf %s %s %s", &wall, &cpu, host,
						 usageList->usage[i].version,
						 usageList->usage[i].threads);
			strncpy(usageList->usage[i].description,&(buff[79]),128);

			usageList->usage[i].wall_total += wall;
			usageList->usage[i].cpu_total  += cpu;
			usageList->usage[i].wall_time   = wall;
			usageList->usage[i].cpu_time    = cpu;
		}
	}

	return usageList;
}


#ifdef TEST_CODE
int main(int argc, char **argv)
{
	long usage;
	char *topdir;

	if (argc != 2) {
		fprintf(stderr,"Usage: du <dirname>\n");
		return -1;
	}
	topdir = argv[1];

	/* Get the disk usage in 512 byte blocks */
	usage = topGetUsage(topdir);

	printf("%lld\n",usage);
}
#endif
