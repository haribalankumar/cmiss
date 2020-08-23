/*
 *  Resource.h
 *
 *    Module for calculating user resources
 */


#ifndef RESOURCE_H
#define RESOURCE_H

#define _POSIX_SOURCE

/* -- Module Datatypes --------------------------------------------------------
*/

typedef struct
{
	double  wall_total;
	double  cpu_total;
	double  wall_time;
	double  cpu_time;
	char    hostname[9];
	char    start_time[25];
	char    version[11];
	char    threads[9];
	char    description[128];
}
USAGE;

typedef struct
{
  int nhost;
	USAGE *usage;
}
USAGELIST;


/* -- Public prototypes --------------------------------------------------------
*/

long long getDiskUsage(char *nodeName);
void fprintDiskUsage(FILE *file, long long diskUsage);

int getCPUUsageByHost(char *filename, char *hostname, double *wall_usage, double *cpu_usage);

USAGELIST *allocUsageList(int nhost);
void freeUsageList(USAGELIST *usageList);
USAGELIST *getUsageFromLog(char *filename, HOSTLIST *hostList);


#endif
