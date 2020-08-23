/*
 *  Directory.c
 *
 *    Module for listing files and directories
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: printf() */

#include <stdlib.h>
    /* For: NULL, malloc(), free() */

#include <string.h>
    /* For: strlen() */

#include <time.h>
    /* For: time_t */

#include <dirent.h>
    /* For: DIR, opendir(), readdir(), etc [System V Flavoured] */

#include <sys/stat.h>
    /* For: stat() */

#include <errno.h>
    /* For: error stuff */

#include "print-error.h"
    /* For: printErrorMessage() */

#include "directory.h"
    /* For: Our public interface */



/* -- Module Datatypes --------------------------------------------------------
*/

#if !defined( BOOLEAN )
#define BOOLEAN int
#define True 1
#define False 0
#endif



/* -- Module Constants --------------------------------------------------------
*/ 

/* 
 *  Max path length, and max file size - different
 *  in different flavours of UNIX. We'll use our
 *  own constants here, but base them on whats available.
 *  POSIX uses (I think) PATH_MAX and FILE_MAX, but
 *  here we'll just include <sys/param.h> and leverage
 *  off MAXPATHLEN and MAXNAMELEN.
 */

/* #include <sys/param.h> */
/* const size_t MaxPathLength     = MAXPATHLEN; */
/* const size_t MaxFilenameLength = MAXNAMELEN; */
#define MaxPathLength 1024
#define MaxFilenameLength 256

extern int errno;


/* -- Method Prototypes -------------------------------------------------------
*/

void addFileInfoToDirectory( FILEINFO *fileInfo, DIRECTORY *directory );



/* -- Test Code ---------------------------------------------------------------
*/ 

void test( char *pathName )
  {   
  DIR           *directory;
  struct dirent *file;

  directory = opendir( pathName );
  if( directory == NULL )
    {
    printf( "some error message here" );
    return;
    }

  file = readdir( directory );
  while( file != NULL )
    {
    printf( "file = %s\n", file->d_name );
    
    file = readdir( directory );
    }

  closedir( directory );
  }


/* -- Private Methods ---------------------------------------------------------
*/

/*
 *  dup()
 */
static char *newStringCopy( const char *source )
  {
  char *tempString;

  tempString = malloc( strlen( source ) + 1 );
  if( tempString == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newStringCopy()\n" );
    return NULL;
    }

  strcpy( tempString, source );
  
  return tempString;
  }


BOOLEAN stringEndsWith( const char *string, const char *end )
  {
  int stringLength;
  int endLength;

  stringLength = strlen( string );
  endLength = strlen( end );

  if( endLength > stringLength ) 
    return False;

  if( strcmp( string + stringLength - endLength, end ) == 0 )
    return True;

  return False;
  }


#define NoDescription      " "
#define ComDescription     "This is a CMISS command file"
#define BasisDescription   "This is a CMISS basis definition file"
#define ElementDescription "This is a CMISS element definition file"
#define EquationDescription "This is a CMISS equation specification file"
#define InitDescription     "This is the initialisation file for CMISS"
#define MaterialDescription "This is a CMISS material definition file"
#define NodeDescription     "This is a CMISS file containing nodal values"
#define SolveDescription    "This is the CMISS algorithm control file"
#define AnalyticDescription "This is a CMISS analytic solutions file"
#define FibreDescription    "This is a CMISS fibre description file"
#define SimpleForm          "This is a simplified parameter entry form"

char *newDescriptionString( const char *path )
  {
  if( stringEndsWith( path, "com" ) == True )
    return newStringCopy( ComDescription );

  if( stringEndsWith( path, "ipbase" ) == True )
    return newStringCopy( BasisDescription );

  if( stringEndsWith( path, "ipelem" ) == True )
    return newStringCopy( ElementDescription );

  if( stringEndsWith( path, "ipequa" ) == True )
    return newStringCopy( EquationDescription );
  
  if( stringEndsWith( path, "ipinit" ) == True )
    return newStringCopy( InitDescription );
  
  if( stringEndsWith( path, "ipmate" ) == True )
    return newStringCopy( MaterialDescription );
 
  if( stringEndsWith( path, "ipnode" ) == True )
    return newStringCopy( NodeDescription );

  if( stringEndsWith( path, "ipsolv" ) == True )
    return newStringCopy( SolveDescription );

  if( stringEndsWith( path, "ipanal" ) == True )
    return newStringCopy( AnalyticDescription );

  if( stringEndsWith( path, "ipfibr" ) == True )
    return newStringCopy( FibreDescription );

  if( stringEndsWith( path, "form" ) == True )
    return newStringCopy( SimpleForm );

  return newStringCopy( NoDescription );
  }


/*
 *  System Specific (i.e. this only work under UNIX) 
 */
static char *newLeafString( const char *path )
  {
  char *tempString;
  char *ptr;

  /* Step 1. find leaf part of name */
  ptr = strrchr( path, '/'  );
  if( ptr == NULL )
    ptr = (char *) path;  /* path in const, ptr is not */
  else
    ptr += 1; /* skip the path element */

  /* Step 2. make a copy for our caller */
  tempString = malloc( strlen( ptr ) + 1 );
  if( tempString == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newLeafString()\n" );
    return NULL;
    }

  strcpy( tempString, ptr );
  
  return tempString;
  }



/* -- Constructors and Destructors --------------------------------------------
*/

FILEINFO *newFileInfo( const char *path )
  {
  FILEINFO    *tempFileInfo;
  struct stat  info;
  int          errorStatus;

  errorStatus = stat( path, &info );

  if( errorStatus == -1 )
    {
    fprintf( stderr, "stat() failed in newFileInfo()\n" );
    return NULL;
    }

  tempFileInfo = malloc( sizeof( FILEINFO ) );
  if( tempFileInfo == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newFileInfo()\n" );
    return NULL;
    }

  tempFileInfo->name = newLeafString( path );
  if( tempFileInfo->name == NULL )
    {
    fprintf( stderr, "newLeafString() returned NULL in newFileInfo()\n" );
    free( tempFileInfo );
    return NULL;
    }

  tempFileInfo->size = info.st_size;

  /* Modification Time */
  tempFileInfo->modificationDate = info.st_mtime;
  
  if( info.st_mode & S_IFREG ) /* Is this a regular file */
    {
    /*
     *  Note: no longer raw octal -- may as well use <sys/stat.h> types, even
		 *  if some comp sci types don't understand them.
     */
    switch( info.st_mode & S_IRWXU )
      {
      case S_IRUSR|S_IWUSR:
        tempFileInfo->access = ReadWriteAccess;
        break;

      case S_IRUSR:
        tempFileInfo->access = ReadOnlyAccess;
        break;

      case S_IWUSR:
        tempFileInfo->access = WriteOnlyAccess;
        break;

      case 00000:
        tempFileInfo->access = NoAccess;
        break;
      }
    }
  else
    {
    tempFileInfo->access = NoAccess;
    }

  if( info.st_mode & S_IFDIR ) /* Is this a directory */
    tempFileInfo->type = DirectoryType;
  else
    tempFileInfo->type = NormalFileType; /* well, close enough for now */

  tempFileInfo->description = newDescriptionString( path );
  if( tempFileInfo->description == NULL )
    {
    fprintf( stderr, "newDescriptionString() returned NULL in "
      "newFileInfo()\n" );
    free( tempFileInfo );
    return NULL;
    }
  
  return tempFileInfo;
  }


void disposeFileInfo( FILEINFO *fileInfo )
  {
  if( fileInfo->name != NULL )
    free( fileInfo->name );

  if( fileInfo->description != NULL )
    free( fileInfo->description );

  free( fileInfo );
  }


FILEINFOITEM *newFileInfoItem( FILEINFO *fileInfo )
  {
  FILEINFOITEM *tempFileInfoItem;

  tempFileInfoItem = malloc( sizeof( FILEINFOLIST ) );
  if( tempFileInfoItem == NULL )
    {
    fprintf( stderr, "Memoey allocation failure in newFileInfoItem()\n" );
    return NULL;
    }

  tempFileInfoItem->fileInfo = fileInfo;
  tempFileInfoItem->nextFileInfoItem = NULL;

  return tempFileInfoItem;
  }


void disposeFileInfoItem( FILEINFOITEM *fileInfoItem )
  {
  disposeFileInfo( fileInfoItem->fileInfo );
  free( fileInfoItem );
  }


void disposeFileInfoItems( FILEINFOITEM *fileInfoItem )
  {
  if( fileInfoItem != NULL )
    {
    disposeFileInfoItems( fileInfoItem->nextFileInfoItem );
    disposeFileInfoItem( fileInfoItem );
    }
  }


FILEINFOLIST *newFileInfoList( void )
  {
  FILEINFOLIST *tempFileInfoList;

  tempFileInfoList = malloc( sizeof( FILEINFOLIST ) );
  if( tempFileInfoList == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newFileInfoList()\n" );
    return NULL;
    }

  tempFileInfoList->firstFileInfoItem = NULL;

  return tempFileInfoList;
  }


void disposeFileInfoList( FILEINFOLIST *fileInfoList )
  {
  if( fileInfoList->firstFileInfoItem != NULL )
    disposeFileInfoItems( fileInfoList->firstFileInfoItem );

  free( fileInfoList );
  }


DIRECTORY *newDirectory( const char *path )
  {
  DIRECTORY     *tempDirectory;
  char           fullFileName[ MaxPathLength + MaxFilenameLength + 1 ];
  DIR           *directory;
  struct dirent *file;
  FILEINFO      *fileInfo;

  tempDirectory = malloc( sizeof( DIRECTORY ) );
  if( tempDirectory == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newDirectory()\n" );
    return NULL;
    }

  tempDirectory->path = newStringCopy( path );  
  if( tempDirectory->path == NULL )
    {
    fprintf( stderr, "newStringCopy() returned NULL in newDirectory()\n" );
    free( tempDirectory );
    return NULL;
    }

  if( path == NULL )
    {
    fprintf( stderr, "path is NULL in newDirectory()\n" );
    free( tempDirectory->path );
    free( tempDirectory );
    return NULL;
    }
  
  directory = opendir( path );
  if( directory == NULL )
    {
		switch (errno)
			{
			case ENOENT:
				{
				printErrorMessage( "newDirectory: Error", "No such directory" );
				exit(0);
				}
			case ENOTDIR:
				{
				printErrorMessage( "newDirectory: Error", "Entry not a directory" );
				exit(0);
				}
			case EACCES:
				{
				printErrorMessage( "newDirectory: Error", "No read permission on directory" );
				exit(0);
				}
			default:
				{
				printErrorMessage( "newDirectory: Error", strerror(errno) );
				exit(0);
				}
			}
    free( tempDirectory->path );
    free( tempDirectory );
    return NULL;
    }

  tempDirectory->fileInfoList = newFileInfoList();
  if( tempDirectory->fileInfoList == NULL )
    {
    fprintf( stderr, "newFileInfoList() returned NULL in newDirectory()\n" );
    free( tempDirectory->path );
    free( tempDirectory );
    return NULL;
    }

  file = readdir( directory );
  while( file != NULL )
    {
    sprintf( fullFileName, "%s/%s", path, file->d_name );

    fileInfo = newFileInfo( fullFileName );
    if( fileInfo != NULL )
      {
      addFileInfoToDirectory( fileInfo, tempDirectory );
      }
    else
      {
      fprintf( stderr, "newFileInfo( %s ) failed in newDirectory() "
        " - skipping\n", fullFileName );
      }
    
    file = readdir( directory );
    }

  closedir( directory );

  return tempDirectory;
  }


void disposeDirectory( DIRECTORY *directory )
  {
  if( directory->path != NULL )   
    free( directory->path );

  if( directory->fileInfoList != NULL )
    disposeFileInfoList( directory->fileInfoList );

  free( directory );
  }



/* -- Methods -----------------------------------------------------------------
*/

void addFileInfoToDirectory( FILEINFO *fileInfo, DIRECTORY *directory )
  {
  FILEINFOITEM *tempFileInfoItem;

  tempFileInfoItem = newFileInfoItem( fileInfo );
  if( tempFileInfoItem == NULL )
    {
    fprintf( stderr, "newFileInfoItem() returned NULL in "
      "addFileInfoToDirectory()\n" );
    return;
    }
  
  /* link into list */
  tempFileInfoItem->nextFileInfoItem =
    directory->fileInfoList->firstFileInfoItem;
  directory->fileInfoList->firstFileInfoItem = tempFileInfoItem;
  }



/* -- Test Code ---------------------------------------------------------------
*/

void printDirectory( DIRECTORY *directory )
  {
  FILEINFO     *fileInfo;
  FILEINFOITEM *fileInfoItem;
printf( "<P>printDirectory()<P>" ); fflush( stdout );

if( directory->fileInfoList == NULL )
{
fprintf( stderr, "directory->fileInfoList == NULL in printDirectory\n" );
return;
}

  fileInfoItem = directory->fileInfoList->firstFileInfoItem;
  while( fileInfoItem != NULL )
    {
    fileInfo = fileInfoItem->fileInfo;

    printf( "%s<BR>", fileInfo->name );

    fileInfoItem = fileInfoItem->nextFileInfoItem;
    }
  }

