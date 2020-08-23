/*
 *  CGI-Decode.c
 *
 *    A module for decoding CGI informaion.
 */


/* -- Standard includes and defines -------------------------------------------
*/

#include <stdlib.h>
    /* For: getenv(), atod(), atoi() */

#include <string.h>
    /* For: strlen() */

#include <ctype.h>
    /* For: isspace() */

#include <stdio.h>
    /* For: printf() */

#include "cgi-decode.h"
    /* For: our public interface */



/* -- Private Data Structures -------------------------------------------------
*/

struct CGIENTRY   /* Completed structure for incomplete type in header */
  {
  char     *name;
  char     *value;
  CGIENTRY *nextEntry;
  };


typedef struct DATAITEM
  {
  int              bufferSize;
  char            *buffer;
  struct DATAITEM *nextItem;
  }
DATAITEM;


typedef struct
  {
  long      dataSize;
  DATAITEM *firstItem;
  DATAITEM *lastItem;
  }
DATALIST;



/* -- Constants and Macros ----------------------------------------------------
*/

/* #define DefaultDataSize 65536 */
#define DefaultDataSize 99999

#define NoAction 0
#define PlusToSpace 1
#define UnescapeChar 2



/* -- Module wide data --------------------------------------------------------
*/

static int actionTable[256];
static int hexTable[256];



/* -- Private method prototypes -----------------------------------------------
*/

void reduceString( char *reducedString, const char *string );
BOOLEAN fuzzyMatch( const char *first, const char *second );

void insertEntry( CGI *cgi, char *name, char *value );

void extractArgs( CGI *cgi, int argc, char *argv[] );
void extractPostFields( CGI *cgi );
void extractGetFields( CGI *cgi );

void copyDataToString( DATALIST *theData, char *theString );
void extendDataList( DATALIST *theList, int extendAmount );
void getStdinData( DATALIST *theDataList );



/* -- Constructors and Destructors --------------------------------------------
*/

CGIENTRY *newCGIEntry( char *name, char *value )
  {
  CGIENTRY *tempEntry;

  tempEntry = malloc( sizeof( CGIENTRY ) );
  if( tempEntry == NULL )
    {
    /* An error would be nice, but where could we print it? */
    return NULL;
    }

  tempEntry->nextEntry = NULL;

  tempEntry->name = name;
  tempEntry->value = value;

  return tempEntry;
  }


CGI *newCGI( void )
  {
  CGI *tempCGI;

  tempCGI = malloc( sizeof( CGI ) );
  if( tempCGI == NULL )
    {
    /* Error! */
    return NULL;
    }

  tempCGI->queryType = Unknown;

  tempCGI->cgiVersion = 0.0;
  tempCGI->httpVersion = 0.0;

  tempCGI->contentType = NULL;
  tempCGI->contentLength = 0;

  tempCGI->argc = 0;
  tempCGI->argv = NULL;

  /* Internal Details */
  tempCGI->firstEntry = NULL;
  tempCGI->queryString = NULL;
  tempCGI->returnString = NULL;

  return tempCGI;
  }


void disposeCGI( CGI *cgi )
  {
  /* !!!!!!!!!!!!!! Finish this */
  }



/* -- Private Constructors/Destructors ----------------------------------------
*/

DATAITEM *newDataItem( int bufferSize )
  {
  DATAITEM *tempDataItem;

  tempDataItem = malloc( sizeof( DATAITEM ) );
  if( tempDataItem == NULL )
    {
    fprintf( stderr, "Out of memory in newDataItem() [1]\n" );

    return NULL;
    }

  tempDataItem->buffer = malloc( sizeof( char ) * bufferSize );
  if( tempDataItem->buffer == NULL )
    {
    fprintf( stderr, "Out of memory in newDataItem() [2]\n" );

    free( tempDataItem );
    return NULL;
    }

  tempDataItem->bufferSize = bufferSize;
  tempDataItem->nextItem = NULL;

  return tempDataItem;
  }


void disposeDataItem( DATAITEM *theDataItem )
  {
  free( theDataItem->buffer );
  free( theDataItem );
  }


void disposeDataItemList( DATAITEM *theDataItem )
  {
  if( theDataItem != NULL )
    {
    disposeDataItemList( theDataItem->nextItem );
    disposeDataItem( theDataItem );
    }
  }


DATALIST *newDataList( void )
  {
  DATALIST *tempDataList;

  tempDataList = malloc( sizeof( DATALIST ) );
  if( tempDataList == NULL )
    {
    fprintf( stderr, "Out of memory in newDataList() [1]\n" );
    
    return NULL;
    }

  tempDataList->firstItem = NULL;
  tempDataList->lastItem = NULL;
  tempDataList->dataSize = 0;

  return tempDataList;
  }


void disposeDataList( DATALIST *theDataList )
  {
  disposeDataItemList( theDataList->firstItem );
  
  free( theDataList );
  }



/* -- String Helper Functions -------------------------------------------------
*/

void reduceString( char *reducedString, const char *string )
  {
  const char *source      = string;
  char       *destination = reducedString;

  while( *source )
    if( ( !isspace( *source ) ) && ( *source != '-' ) && ( *source != '_' ) )
      *destination++ = toupper( *source++ );
    else
      source++;  

  *destination = '\0'; /* terminate the string */
  }


BOOLEAN fuzzyMatch( const char *first, const char *second )
  {
  char *firstReduced;
  char *secondReduced;
  BOOLEAN  result;

  firstReduced =  malloc( strlen( first ) + 1 );
  secondReduced =  malloc( strlen( second ) + 1 );

  reduceString( firstReduced, first );
  reduceString( secondReduced, second );

  if( strcmp( firstReduced, secondReduced ) == 0 )
    result = True;
  else
    result = False;

  free( firstReduced );
  free( secondReduced );
  
  return result;
  }


/*
 *  Portable, fast, hex decoding table
 * 
 *  Please note that I'm not sure about exactly *what* a character constant
 *  such as 'a' returns - ISO C 6.1.3.4: Character Constants states that
 *  it is an int, but doesn't mention the sign. An explicit cast *shouldn't*
 *  be nessessary in the next two function, but who know what sort of
 *  bizarre character stes there are out there...
 */
void initHexTable( int ht[256] )
  {
  int i;

  for( i = 0; i < 256; i++ )
    ht[i] = 0;

  ht['0']=0; ht['1']=1; ht['2']=2; ht['3']=3; ht['4']=4;
  ht['5']=5; ht['6']=6; ht['7']=7; ht['8']=8; ht['9']=9;
  ht['a']=10; ht['b']=11; ht['c']=12; ht['d']=13; ht['e']=14; ht['f']=15;
  ht['A']=10; ht['B']=11; ht['C']=12; ht['D']=13; ht['E']=14; ht['F']=15; 
  }


void initActionTable( int actionTable[256] )
  {
  int i;

  for( i = 0; i < 256; i++ )
    actionTable[i] = NoAction;

  actionTable['+'] = PlusToSpace;
  actionTable['%'] = UnescapeChar;
  }


/*
 *  Decodes an 'Content-Type: x-url-encoded' encoded string in place
 */ 
void decode( char *theString )
  {
  unsigned char *source;
  unsigned char *dest;  /* A *must* for using a char to dereference arrays */
  unsigned char  rawChar;

  source = (unsigned char *) theString;
  dest = (unsigned char *) theString;

  while( *source != '\0' )
    {
    switch( actionTable[*source] )
      {
      case PlusToSpace:
        *dest++ = ' ';
        source++; 
        break;

      case UnescapeChar:
        source++;
        rawChar = hexTable[*source++] * 16;
        rawChar += hexTable[*source++];
        *dest++ = rawChar;
        break;

      default:
        *dest++ = *source++;
      }
    }

  *dest = '\0'; /* Terminate our destination string */
  }

char *newStringFromBounds( char *start, char *end )
  {
  char *tempString;
  int   length;
  int   i;

  length = end - start;

  tempString = malloc( sizeof( char ) * (length + 1) );
  if( tempString == NULL )
    {
    fprintf( stderr, "memory allocation failure in newStringFromBounds()\n" );
    return NULL;
    }

  for( i = 0; i < length; i++ )
    tempString[i] = start[i];

  tempString[length] = '\0';

  return tempString;
  }



/* -- External CGI Library Calls ----------------------------------------------
*/

BOOLEAN nameExists( CGI *cgi, char *name )
  {
  CGIENTRY *entry;

  entry = cgi->firstEntry;

  while( entry != NULL )
    {
    if( fuzzyMatch( entry->name, name ) )
      return True;
    
    entry = entry->nextEntry;
    }

  return False;
  }


char *lookupString( CGI *cgi, char *name )
  {
  CGIENTRY *entry;

  entry = cgi->firstEntry;

  while( entry != NULL )
    {
    if( fuzzyMatch( entry->name, name ) )
      return entry->value;

    entry = entry->nextEntry;
    }

  /* Perhaps this should give an error instead? */
  return NULL;
  }


int lookupNumber( CGI *cgi, char *name )
  {
  CGIENTRY *entry;

  entry = cgi->firstEntry;

  while( entry != NULL )
    {
    if( fuzzyMatch( entry->name, name ) )
      return atoi( entry->value );

    entry = entry->nextEntry;
    }

  /* Perhaps this should give an error instead? */
  return 0;
  }



/* -- Private Methods ---------------------------------------------------------
*/

void copyDataToString( DATALIST *theData, char *theString )
  {
  char     *currentPtr;
  DATAITEM *currentItem;
  int       amountLeft;

  currentPtr = theString;
  currentItem = theData->firstItem;

  amountLeft = theData->dataSize;

  while( currentItem != NULL )
    {
    if( amountLeft < currentItem->bufferSize )
      {
      memcpy( currentPtr, currentItem->buffer, amountLeft );
      break;
      }
    else
      {
      memcpy( currentPtr, currentItem->buffer, currentItem->bufferSize );
      amountLeft -= currentItem->bufferSize;
      currentPtr += currentItem->bufferSize;
      }

    currentItem = currentItem->nextItem;
    }
  
  theString[theData->dataSize] = '\0'; /* Zero terminate string for caller */
  }


void extendDataList( DATALIST *theList, int extendAmount )
  {
  DATAITEM *tempDataItem;

  tempDataItem = newDataItem( extendAmount );
  if( tempDataItem == NULL )
    {
    fprintf( stderr, "newDataItem() failed in extendDataList()\n" );
    return;
    }

  if( theList->lastItem == NULL )
    {
    theList->firstItem = tempDataItem;
    theList->lastItem = tempDataItem;
    }
  else
    {
    theList->lastItem->nextItem = tempDataItem;
    theList->lastItem = tempDataItem;
    }
  }


void getStdinData( DATALIST *theDataList )
  {
  int       contentLength;
  char     *lengthString;
  int       amountRead;   /* NB sizeof(char) == 1  [ISO C 6.3.3.4] */ 

  lengthString = getenv ( "CONTENT_LENGTH" );

  if( lengthString != NULL )
    {
    /* lengthString contains the ascii representation of the Content-Length: */

    contentLength = atoi( lengthString );

    extendDataList( theDataList, contentLength );

    amountRead = fread( theDataList->lastItem->buffer, sizeof( char ),
                        contentLength, stdin );

    theDataList->dataSize = amountRead;
    }
  else
    {
    /* We don't know *how* much data there is to come */

    extendDataList( theDataList, DefaultDataSize );

    amountRead = fread( theDataList->lastItem->buffer, sizeof( char ),
                        contentLength, stdin );

    theDataList->dataSize += amountRead;

    while( amountRead == DefaultDataSize )
      {
      extendDataList( theDataList, DefaultDataSize );
      
      amountRead = fread( theDataList->lastItem->buffer, sizeof( char ),
                        contentLength, stdin );

      theDataList->dataSize += amountRead;
      }
    }
  }


void insertEntry( CGI *cgi, char *name, char *value )
  {
  CGIENTRY* tempEntry;

  if( nameExists( cgi, name ) == False )
    {
    tempEntry = newCGIEntry( name, value );
    if( tempEntry == NULL )
      {
      /* Error! Maybe give details and exit later? */
      return;
      }
    
    tempEntry->nextEntry = cgi->firstEntry;
    cgi->firstEntry = tempEntry;
    }
  else
    {
    /* Error, but where to report it? */
    return;
    }

  return;  /* Normal return, no errors */
  }


void extractArgs( CGI *cgi, int argc, char *argv[] )
  {
  /*
   *  Hmm, do nothing for now, maybe do this later (do we need argc/argv??) 
   *
  printf( "extractArgs() called\n" );
   *
   */
  }


/*
 *  This function makes my skin crawl. I think it must be time to invest
 *  in a decent RE engine in preference to doing the abonimation below
 *  any more times...
 */
BOOLEAN mapQuery( char *queryString )
  {
  char *ptr;
  int   state;

  ptr = queryString;
  state = 0;

  while( *ptr != '\0' )
    {
    switch( state ) 
      {
      case 0:
        if( isdigit( *ptr ) )
          state = 1;
        else
          return False;
        break;

      case 1:
        if( *ptr == ',' )
          state = 2;
        else
          if( isdigit( *ptr ) == 0 )
            return False;
        break;

      case 2:
        if( isdigit( *ptr ) )
          state = 3;
        else
          return False;
        break;

      case 3:
        if( *ptr == '\0' || *ptr == '\n' || *ptr == ' ' )
          return True;

        if( isdigit( *ptr ) == 0 )
          return False;

        break;
      }
 
    ptr++;
    }

  if( state == 3 )
    return True;
  else
    return False;
  }


/* 
 *  This could be tidier as well, but I guess it will do.
 */
void decodeMapQuery( CGI *cgi )
  {
  char *startPosition;
  char *commaPosition;
  char *endPosition;
  char *name;
  char *value;
  char *x = "x";
  char *y = "y";

  startPosition = cgi->queryString;

  commaPosition = strchr( cgi->queryString, ',' );
  if( commaPosition == NULL )
    {
    fprintf( stderr, "Incorrectly formed input string in decodeMapQuery()\n" );
    return;
    }

  endPosition = cgi->queryString + strlen( cgi->queryString );

  name = newStringFromBounds( x, x + strlen( x ) );
  value = newStringFromBounds( startPosition, commaPosition );

  insertEntry( cgi, name, value );

  name = newStringFromBounds( y, y + strlen( y ) );
  value = newStringFromBounds( commaPosition + 1, endPosition );

  insertEntry( cgi, name, value );
  }


void extractPairs( CGI *cgi )
  {
  char    *name;
  char    *value;
  char    *currentPtr;
  char    *equalsPtr;
  char    *ampersandPtr;
  BOOLEAN morePairsLeft;


  if( strchr( cgi->queryString, '=' ) == NULL )
    {
    if( mapQuery( cgi->queryString ) == True )
      {
      decodeMapQuery( cgi );
      }
    else
      {
      /* There are no name/value pairs - unescape just queryString instead */
      decode( cgi->queryString );
      }
    return;
    }

  /* Otherwise go through extracting pairs name=value sperated by &'s */
  currentPtr = cgi->queryString;
  morePairsLeft = True;
  
  while( morePairsLeft )
    {
    equalsPtr = strchr( currentPtr, '=' );
    ampersandPtr = strchr( currentPtr, '&' );

    if( equalsPtr == NULL )
      {
      /* This is not a pair, the queryString is poorly formed */
      fprintf( stderr, "poorly formed query string found in extractPairs()\n" );
      return;
      }

    if( ampersandPtr == NULL )
      {
      /* No more pairs after this one */
      morePairsLeft = False;
      ampersandPtr = currentPtr + strlen( currentPtr );/* last char + 1 */
      }

    name = newStringFromBounds( currentPtr, equalsPtr );
    value = newStringFromBounds( equalsPtr + 1, ampersandPtr );

    decode( name );
    decode( value );

    /*printf( "name = %s, value = %s\n", name, value );*/
  
    insertEntry( cgi, name, value );

    currentPtr = ampersandPtr + 1;
    }
  }


/*
 *  The POST method delivers input via a stdin pipe, to the program, along
 *  with certain MIME data specified in environment variables (such as
 *  Content-type: and Content-length:)
 */
void extractPostFields( CGI *cgi )
  {
  char     *queryString;
  int       queryLength;
  DATALIST *theData;

  theData = newDataList();
  getStdinData( theData );

  queryLength = theData->dataSize + 1;
  queryString = malloc( sizeof( char ) * queryLength );
  copyDataToString( theData, queryString );
  cgi->queryString = queryString;

  disposeDataList( theData );

  extractPairs( cgi );
  }


/*
 *  The GET method delivers input via the environment variable  
 *  "QUERY_STRING".
 */
void extractGetFields( CGI *cgi )
  {
  char *queryEnvironment;
  char *queryString;
  int   queryLength;

  queryEnvironment = getenv( "QUERY_STRING" );
  if( queryEnvironment == NULL )
    {
    /* Hmmm, no query string, make like it's a zero length one. */
    queryLength = 1;
    queryString = malloc( queryLength );
    strcpy( queryString, "" );
    cgi->queryString = queryString;
    }
  else
    {
    queryLength = strlen( queryEnvironment ) + 1;
    queryString = malloc( queryLength );
    strcpy( queryString, queryEnvironment ); /* Make local copy for ourselves */
    cgi->queryString = queryString;
    }

  extractPairs( cgi );
  }


/*
getNextChar( INPUT *theInput )
  {
  if( theInput( ))
  }
*/


CGI *getCGIEnvironment( int argc, char *argv[] )
  {
  char *cgiString;
  char *httpString;
  char *queryType;
  CGI  *cgi;

  cgi = newCGI();

  initActionTable( actionTable );
  initHexTable( hexTable );

  /*
   *  If there were any command line arguments, extract them
   */
  if( argc > 0 )
    extractArgs( cgi, argc, argv );
  else
    cgi->argc = 0;

  /*
   *  What method was used for this request (GET or POST only)
   */ 
  queryType = getenv( "REQUEST_METHOD" );
  if( queryType == NULL )
    {
    fprintf( stderr, "REQUEST_METHOD not defined in getCGIEnvironment()\n" );

    disposeCGI( cgi );
    return NULL;
    }

  if( fuzzyMatch( queryType, "POST" ) == True )
    {
    cgi->queryType = Post;
    }
  else
    {
    if( fuzzyMatch( queryType, "GET" ) == True )
      {
      cgi->queryType = Get;
      }
    else
      {
      cgi->queryType = Unknown;
      }
    }

  /*
   *  Extract the contents of the query string, if it exists
   */
  switch( cgi->queryType )
    {
    case Post:
      extractPostFields( cgi );
      break;

    case Get:
      extractGetFields( cgi );
      break;

    case Unknown:
      fprintf( stderr, "Unknown query type in getCGIEnvironment()\n" );
      break;
    }

  /*
   *  Extract the version information from the CGI and HTTP protcols
   */
  cgiString = getenv( "GATEWAY_INTERFACE" );
  if( cgiString != NULL )
    {
    char *lastSlash = strrchr( cgiString, '/' );
    
    if( lastSlash != NULL )
      {
      cgi->cgiVersion = atof( lastSlash + 1 );
      }
    }

  httpString = getenv( "SERVER_PROTOCOL" );
  if( httpString != NULL )
    {
    char *lastSlash = strrchr( httpString, '/' );
    
    if( lastSlash != NULL )
      {
      cgi->httpVersion = atof( lastSlash + 1 );
      }
    }

  return cgi;
  }

