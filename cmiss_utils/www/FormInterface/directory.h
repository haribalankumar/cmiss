/*
 *  Directory.h
 *
 *    Public interface for the Directory Module
 */
#ifndef DIRECTORY_H
#define DIRECTORY_H


/* -- Required Include Files --------------------------------------------------
*/

#include <time.h>
    /* For: time_t */



/* -- Module Constants --------------------------------------------------------
*/ 

/* Access */
#define NoAccess        0
#define ReadOnlyAccess  1
#define WriteOnlyAccess 2
#define ReadWriteAccess 3

/* Type */
#define NormalFileType 0
#define DirectoryType  1



/* -- Module Data Structures --------------------------------------------------
*/

typedef struct
  {
  char   *name;
  int     size;
  int     access;
  int     type;
  time_t  modificationDate;
  char   *description;
  }
FILEINFO;


typedef struct FILEINFOITEM
  {
  FILEINFO *fileInfo;
  
  struct FILEINFOITEM *nextFileInfoItem;
  }
FILEINFOITEM;


typedef struct 
  {
  int           numberOfItems;
  FILEINFOITEM *firstFileInfoItem;
  }
FILEINFOLIST;


typedef struct
  {
  FILEINFOLIST *fileInfoList;
  char         *path;
  }
DIRECTORY;



/* -- Constructors and Destructors --------------------------------------------
*/

FILEINFO *newFileInfo( const char *path );
void disposeFileInfo( FILEINFO *fileInfo );

FILEINFOITEM *newFileInfoItem( FILEINFO *fileInfo );
void disposeFileInfoItem( FILEINFOITEM *fileInfoItem );
void disposeFileInfoItems( FILEINFOITEM *fileInfoItem );

DIRECTORY *newDirectory( const char *path );
void disposeDirectory( DIRECTORY *directory );



/* -- Module Public Methods ---------------------------------------------------
*/

void printDirectory( DIRECTORY *directory );
void filterDirectory( DIRECTORY *directory, void *filterPatterns );
void sortDirectory( DIRECTORY *directory );


#endif

