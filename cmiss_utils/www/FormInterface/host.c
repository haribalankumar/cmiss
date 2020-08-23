/*
 *  Host.c
 *
 *    Manages Hosts
 */


/* -- Include Directives ------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stderr */

#include <stdlib.h>
    /* For: malloc(), free(), NULL, srand(), rand() */

#include <string.h>
    /* For: strcpy() */

#include "groups.h"
    /* For: GROUP, GROUPFILE, setGroupFilePos(), etc, etc */

#include "rules.h"
    /* For: RULES, newRules(), etc */

#include "rule-parse.h"
    /* For: parseRules() */

#include "print-error.h"
    /* For: parseErrorMessage() */

#include "host.h"
    /* For: our public methods and datatypes */



/* -- Private Module Constants ------------------------------------------------
*/

/* String Constants */
#define StringMatch 0



/* -- Private Module Utilities ------------------------------------------------
*/

static void stripTrailingSpace( char *string )
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

HOST *newHost( char *hostname, char *username, char *os, char *ncpu, 
	char *homeDirectory )
  {
  HOST *tempHost;

  tempHost = malloc( sizeof( HOST ) );
  if( tempHost == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newHost() [1]\n" );
    return NULL;
    }

  tempHost->hostname = malloc( strlen( hostname ) + 1 );
  if( tempHost->hostname == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newHost() [2]\n" );
    free( tempHost );
    return NULL;
    }

  tempHost->username = malloc( strlen( username ) + 1 );
  if( tempHost->username == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newHost() [3]\n" );
    free( tempHost->hostname );
    free( tempHost );
    return NULL;
    }

  tempHost->os = malloc( strlen( os ) + 1 );
  if( tempHost->os == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newHost() [4]\n" );
    free( tempHost->username );
    free( tempHost->hostname );
    free( tempHost );
    return NULL;
    }

  tempHost->ncpu = malloc( strlen( ncpu ) + 1 );
  if( tempHost->ncpu == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newHost() [5]\n" );
    free( tempHost->os );
    free( tempHost->username );
    free( tempHost->hostname );
    free( tempHost );
    return NULL;
    }

  tempHost->homeDirectory = malloc( strlen( homeDirectory ) + 1 );
  if( tempHost->homeDirectory == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newHost() [6]\n" );
    free( tempHost->ncpu );
    free( tempHost->os );
    free( tempHost->username );
    free( tempHost->hostname );
    free( tempHost );
    return NULL;
    }

  strcpy( tempHost->hostname, hostname );
  strcpy( tempHost->username, username );
  strcpy( tempHost->os, os );
  strcpy( tempHost->ncpu, ncpu );
  strcpy( tempHost->homeDirectory, homeDirectory );

  stripTrailingSpace( tempHost->hostname );
  stripTrailingSpace( tempHost->username );
  stripTrailingSpace( tempHost->os );
  stripTrailingSpace( tempHost->ncpu );
  stripTrailingSpace( tempHost->homeDirectory );

  return tempHost;
  }


HOST *duplicateHost( HOST *host )
  {
  HOST *tempHost;

  tempHost = newHost( host->hostname, host->username, host->os, 
										  host->ncpu, host->homeDirectory );
  if( tempHost == NULL )
    {
    fprintf( stdout, "newHost() returned NULL in duplicateHost()\n" );
    return NULL;
    }

  return tempHost;
  }


void disposeHost( HOST *host )
  {
  if( host->hostname != NULL )
    free( host->hostname );
  if( host->username != NULL )
    free( host->username );
  if( host->os != NULL )
    free( host->os );
  if( host->ncpu != NULL )
    free( host->ncpu );
  if( host->homeDirectory != NULL )
    free( host->homeDirectory );
  free( host );
  }


HOSTITEM *newHostItem( HOST *host )
  {
  HOSTITEM *tempHostItem;
  
  tempHostItem = malloc( sizeof( HOSTITEM ) );
  if( tempHostItem == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newHostItem()\n" );
    return NULL;
    }

  tempHostItem->host = host;
  tempHostItem->nextHostItem = NULL;

  return tempHostItem;
  }


void disposeHostItem( HOSTITEM *hostItem )
  {
  disposeHost( hostItem->host );
  free( hostItem );
  }


void disposeHostItems( HOSTITEM *hostItem )
  {
  if( hostItem != NULL )
    {
    disposeHostItems( hostItem->nextHostItem );
    disposeHostItem( hostItem );
    }
  }


HOSTLIST *newHostList( void )
  {
  HOSTLIST *tempHostList;

  tempHostList = malloc( sizeof( HOSTLIST ) );
  if( tempHostList == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newHostList()\n" );
    return NULL;
    }

  tempHostList->firstHostItem = NULL;

  return tempHostList;
  }


void disposeHostList( HOSTLIST *hostList )
  {
  if( hostList->firstHostItem != NULL )
    disposeHostItems( hostList->firstHostItem );

  free( hostList );
  }



/* -- Public Methods ----------------------------------------------------------
*/

HOST *lookupHost( char *hostname, HOSTLIST *hostList )
  {
  HOSTITEM *hostItem;
  HOST     *host;

  if( hostname == NULL )
    return NULL;

  hostItem = hostList->firstHostItem;
  while( hostItem != NULL )
    {
    host = hostItem->host;

    if( tagCompare( hostname, host->hostname ) == True )
      return host;
    
    hostItem = hostItem->nextHostItem;
    }

  return NULL;
  }


void addHostToHostList( HOST *host, HOSTLIST *hostList )
  {
  HOSTITEM *hostItem, *hostPtr;

  hostItem = newHostItem( host );
  if( hostItem == NULL )
    {
    fprintf( stderr, "newHostItem() returned NULL in addHostToHostList()\n" );
    return;
    }

#ifdef OLD_CODE
	/* This code reverses the file order */
  hostItem->nextHostItem = hostList->firstHostItem;
  hostList->firstHostItem = hostItem;
#else
	/* This code doesn't reverse the file order */
	hostItem->nextHostItem = NULL;

	if (hostList->firstHostItem == NULL)
		{
		hostList->firstHostItem = hostItem;
	  }
	else
		{
		hostPtr = hostList->firstHostItem;
		while( hostPtr->nextHostItem != NULL )
			{
			hostPtr = hostPtr->nextHostItem;
			}
		hostPtr->nextHostItem = hostItem;
	  }
#endif
  }



void readHostList( HOSTLIST *hostList, char *baseName )
  {
  GROUPFILE *file;
  RULES     *rules;
  GROUP     *group;
  TAG       *hostnameTag;
  TAG       *usernameTag;
  TAG       *osTag;
  TAG       *ncpuTag;
  TAG       *homeDirectoryTag;
  HOST      *host;
  char      *string;

	string = malloc( strlen( baseName ) + 7 ); /* + ".rules" + '\0' */
  if( string == NULL )
    {
    fprintf( stderr, "Memory allocation failure in readHostList()\n" );
    return;
    }

  sprintf( string, "%s.list", baseName );
  file = newGroupFile( string );
  if( file == NULL )
    {
    printErrorMessage( "readHostList: Error",
      "newGroupFile() returned NULL in readHostList()" );
    free( string );
    return;
    }
  
  rules = newRules();
  if( rules == NULL )
    {
    printErrorMessage( "readHostList: Error", 
			"newRules() returned NULL in readHostList()" );
    free( string );
    disposeGroupFile( file );
    return;
    }

  sprintf( string, "%s.rules", baseName );
  parseRuleFile( rules, string );
  
  /* -- No error checking performed */

  group = getGroup( file, rules );
  while( group != NULL )
    {
    hostnameTag = lookupTag( group, "hostname" );
    usernameTag = lookupTag( group, "username" );
    osTag = lookupTag( group, "os" );
    ncpuTag = lookupTag( group, "ncpu" );
    homeDirectoryTag = lookupTag( group, "home-directory" );

    if( hostnameTag != NULL && usernameTag != NULL && osTag != NULL && 
			  ncpuTag != NULL && homeDirectoryTag != NULL )
      {
      host = newHost( hostnameTag->body, usernameTag->body, osTag->body, 
        ncpuTag->body, homeDirectoryTag->body );
      addHostToHostList( host, hostList );
      }

    discardGroup( file );
    group = getGroup( file, rules );
    }

  disposeGroupFile( file );
  disposeRules( rules );
  free( string );
  }


void printHostList( FILE *hostFile, HOSTLIST *hostList )
  {
	HOSTITEM *hostItem;
	FILE *file;

	file = hostFile;

	if (hostFile == NULL)
		{
		file = stdout;
		}

	hostItem = hostList -> firstHostItem;
  while( hostItem != NULL )
    {
		fprintf( file, "hostname: %s\n",       hostItem -> host -> hostname );
		fprintf( file, "username: %s\n",       hostItem -> host -> username );
		fprintf( file, "os: %s\n",             hostItem -> host -> os );
		fprintf( file, "ncpu: %s\n",           hostItem -> host -> ncpu );
		fprintf( file, "home-directory: %s\n", hostItem -> host -> homeDirectory );
		hostItem = hostItem -> nextHostItem;
		if (hostItem != NULL) fprintf( file, "\n");
    }

  }



