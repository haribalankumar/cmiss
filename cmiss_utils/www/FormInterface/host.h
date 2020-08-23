/*
 *  Host.h
 *
 *    Public Interface to the Host routines
 */
#ifndef HOST_H
#define HOST_H


/* -- Required Include Files --------------------------------------------------
*/

#include "rules.h"
    /* For: RULES, newRules(), etc */

#include "rule-parse.h"
    /* For: parseRules() */



/* -- Module Datatypes --------------------------------------------------------
*/

typedef struct
  {
  char *hostname;
  char *username;
  char *os;
  char *ncpu;
  char *homeDirectory;
  }
HOST;


typedef struct HOSTITEM
  {
  HOST *host;

  struct HOSTITEM *nextHostItem;
  }
HOSTITEM;


typedef struct 
  {
  HOSTITEM *firstHostItem;
  }
HOSTLIST;



/* -- Module Constructor and Destructors --------------------------------------
*/

HOST *newHost( char *hostname, char *username, char *os, char *ncpu,
	char *homeDirectory );
HOST *duplicateHost( HOST *host );
void disposeHost( HOST *host );

HOSTITEM *newHostItem( HOST *host );
void disposeHostItem( HOSTITEM *hostItem );
void disposeHostItems( HOSTITEM *hostItem );

HOSTLIST *newHostList( void );
void disposeHostList( HOSTLIST *hostList );


/* -- Public Methods ----------------------------------------------------------
*/

HOST *lookupHost( char *hostname, HOSTLIST *hostList );
void addHostToHostList( HOST *host, HOSTLIST *hostList );
void readHostList( HOSTLIST *hostList, char *baseName );
void printHostList( FILE *hostFile, HOSTLIST *hostList );

#endif
