/*
 *  Call Tree.c
 *  
 *    Generates the calling tree information from the CMISS source.
 *
 * Modified by:
 *   Karl Tomlinson 4Nov99:  function call handling.
 *   Karl Tomlinson 8Dec99:  external symbol handling.
 *   Karl Tomlinson 9Oct00:  modified calls exclusion handling and added
 *                   groups for routines that are not documented in CMISS.
 *
 *    "Now a life of leisure, and a pirate's treasure,
 *     Don't make much for tragedy.
 *     But it's a sad man my friend who's livin' in his own skin,
 *     And can't stand the company.
 *     Every fool's got a reason to feelin' sorry for himself,
 *     And turn his heart to stone.
 *     Tonight this fool's halfway to heaven and just a mile outta hell,
 *     And I feel like I'm comin' home.
 *     
 *     These are better days baby,
 *     There's better days shining through,
 *     These are better days,
 *     Better days with a girl like you."
 *                                                  -- Bruce Springsteen
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: FILE, fprintf(), etc */

#include <errno.h>

#include <stdlib.h>
    /* For: NULL, malloc(), free() */

#include <string.h>
    /* For: strlen(), strcat(), strcpy() */

#include <ctype.h>
    /* toupper() */

#include "simple-parser.h"
    /* For: PARSESTATE, newParseState(), skip(), extract(), etc */

#include "lineio.h"
    /* For: FILESTATE, newFileState(), getLine(), discardLine() */

#include "tags.h"
    /* For: TAG, newTag(), etc */

#include "wwwpaths.h"
    /* For CMISS_WWW_ROOT LOOKUP_URLPATH */

#include "groups.h"
    /* For: newGroupFromHeadTag */

#include "collection.h"
    /* For: COLLECTION, newCollectionFromFile(), lookupGroupList(), etc */

#include "reducers.h"
    /* For: reduce*() */



/* -- Module Datatypes --------------------------------------------------------
*/

typedef struct
  {
  FILESTATE *file;
  char      *currentLine;
  BOOLEAN    readPending;
  int        currentLineNumber;  /* conceptual name clash with currentLine */
  int        currentByteOffset;
  } 
FORTRANFILE;


typedef struct STRINGITEM
  {
  char *string;
  
  struct STRINGITEM *nextStringItem;
  }
STRINGITEM;


typedef struct
  {
  STRINGITEM *firstStringItem;
  STRINGITEM *lastStringItem;
  }
STRINGLIST;


typedef struct
  {
  char *name;
  }
CALL;


typedef struct CALLITEM
  {
  CALL *call;

  struct CALLITEM *nextCallItem;
  }
CALLITEM;


typedef struct
  {
  CALLITEM *firstCallItem;
/*    BOOLEAN   truncated;     */
  }
CALLLIST;


typedef struct
  {
  char     *routineName;

  char     *moduleName;
  int       routineStartLine;
  int       routineEndLine;
  int       routineStartPosition;
  int       routineEndPosition;

  CALLLIST *calls;
  CALLLIST *calledFrom;
  }
ROUTINE;


typedef struct ROUTINEITEM
  {
  ROUTINE *routine;

  struct ROUTINEITEM *nextRoutineItem;
  }
ROUTINEITEM;


typedef struct
  {
  ROUTINEITEM *firstRoutineItem;
  }
ROUTINELIST;


typedef struct
  {
  char        *moduleName;
  ROUTINELIST *routineList;
  }
MODULE;


typedef struct MODULEITEM
  {
  MODULE *module;

  struct MODULEITEM *nextModuleItem;
  }
MODULEITEM;


typedef struct
  {  
  MODULEITEM *firstModuleItem;
  }
MODULELIST;



/* -- Module Constants --------------------------------------------------------
*/

/* Files */
char *DefaultSubroutineFileName = CMISS_WWW_ROOT LOOKUP_URLPATH "/subroutines";
char *DefaultFunctionFileName = CMISS_WWW_ROOT LOOKUP_URLPATH "/functions";
char *DefaultModuleBaseName = CMISS_WWW_ROOT LOOKUP_URLPATH "/modules";
char *DefaultModuleListFileName = "../MasterLists/module-list";

char *MasterModuleListFile = CMISS_WWW_ROOT CM_URLPATH
  "/modules-list.html";

/* Excluded routine names from calls tag */
static char *callsExclusionList[] = { "ENTERS", "EXITS", NULL };
/*  static const int maxCalledFromSize = 100; */

/* State Names */
#define InsideSubroutine 1
#define InsideFunction 2
#define OutsideRoutine 4

/* Programmer Amelioration */
#define StringMatch 0



/* -- Constructors and Destructors --------------------------------------------
*/

CALL *newCall( char *name )
  {
  CALL *tempCall;

  tempCall = malloc( sizeof( CALL ) );
  if( tempCall == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newCall() [1]\n" );
    return NULL;
    }

  tempCall->name = malloc( strlen( name ) + 1 );
  if( tempCall->name == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newCall() [2]\n" );
    free( tempCall );
    return NULL;
    }

  strcpy( tempCall->name, name );

  return tempCall;
  }


void disposeCall( CALL *call )
  {
  if( call->name != NULL )
    free( call->name );

  free( call );
  }


CALLITEM *newCallItem( CALL *call )
  {
  CALLITEM *tempCallItem;

  tempCallItem = malloc( sizeof( CALLITEM ) );
  if( tempCallItem == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newCallItem()\n" );
    return NULL;
    }

  tempCallItem->call = call;

  return tempCallItem;
  }


void disposeCallItem( CALLITEM *callItem )
  {
  disposeCall( callItem->call );
  free( callItem );
  }


void disposeCallItems( CALLITEM *callItem )
  {
  if( callItem != NULL )
    {
    disposeCallItems( callItem->nextCallItem );
    disposeCallItem( callItem );
    } 
  }


CALLLIST *newCallList( void )
  {
  CALLLIST *tempCallList;

  tempCallList = malloc( sizeof( CALLLIST ) ); 

  if( tempCallList == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newCallList()\n" );
    return NULL;
    }

  tempCallList->firstCallItem = NULL;
/*    tempCallList->truncated = False; */
  
  return tempCallList;
  }


void disposeCallList( CALLLIST *callList )
  {
  disposeCallItems( callList->firstCallItem );
  free( callList );
  }


ROUTINE *newRoutine( char *name )
  {
  ROUTINE *tempRoutine;

  tempRoutine = malloc( sizeof( ROUTINE ) );
  if( tempRoutine == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newRoutine() [1]\n" );
    return NULL;
    }

  tempRoutine->routineName = malloc( strlen( name ) + 1 );
  if( tempRoutine->routineName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newRoutine() [2]\n" );
    free( tempRoutine );
    }

  tempRoutine->calls = newCallList();
  if( tempRoutine->calls == NULL )
    {
    fprintf( stderr, "newCallList() returned NULL in newRoutine() [1]\n" );
    free( tempRoutine->routineName );
    free( tempRoutine );
    }

  tempRoutine->calledFrom = newCallList();
  if( tempRoutine->calledFrom == NULL )
    {
    fprintf( stderr, "newCallList() returned NULL in newRoutine() [1]\n" );
    disposeCallList( tempRoutine->calls );
    free( tempRoutine->routineName );
    free( tempRoutine );
    }

  strcpy( tempRoutine->routineName, name );
  tempRoutine->moduleName = NULL;
  tempRoutine->routineStartLine = -1;
  tempRoutine->routineEndLine = -1;
  tempRoutine->routineStartPosition = -1;
  tempRoutine->routineEndPosition = -1;

  return tempRoutine;
  }


void disposeRoutine( ROUTINE *routine )
  {
  printf( "Please write disposeRoutine()\n" );
  }


ROUTINEITEM *newRoutineItem( ROUTINE *routine )
  {
  ROUTINEITEM *tempRoutineItem;

  tempRoutineItem = malloc( sizeof( ROUTINEITEM ) );
  if( tempRoutineItem == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newRoutineItem()\n" );
    return NULL;
    }

  tempRoutineItem->routine = routine;
  tempRoutineItem->nextRoutineItem = NULL;

  return tempRoutineItem;
  }


void disposeRoutineItem( ROUTINEITEM *routineItem )
  {
  if( routineItem->routine != NULL )
    disposeRoutine( routineItem->routine );
  free( routineItem );
  }


void disposeRoutineItems( ROUTINEITEM *routineItem )
  {
  if( routineItem != NULL )
    {
    disposeRoutineItems( routineItem->nextRoutineItem );
    disposeRoutineItem( routineItem );
    }
  }


ROUTINELIST *newRoutineList( void )
  {
  ROUTINELIST *tempRoutineList;

  tempRoutineList = malloc( sizeof( ROUTINELIST ) );
  if( tempRoutineList == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newRoutineList()\n" );
    return NULL;
    }

  tempRoutineList->firstRoutineItem = NULL;

  return tempRoutineList;
  }

void disposeRoutineList( ROUTINELIST *routineList )
  {
  printf( "Please write disposeRoutineList()\n" );
  }


MODULE *newModule( char *moduleName )
  {
  MODULE *tempModule;

  tempModule = malloc( sizeof( MODULE ) );
  if( tempModule == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newModule() [1]\n" );
    return NULL;
    }

  tempModule->moduleName = malloc( strlen( moduleName ) + 1 );
  if( tempModule->moduleName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newModule() [2]\n" );
    free( tempModule );
    return NULL;
    }

  tempModule->routineList = newRoutineList();
  if( tempModule->routineList == NULL )
    {
    fprintf( stderr, "newRoutineList() returned NULL in newModule()\n" );
    free( tempModule->moduleName );
    free( tempModule );
    return NULL;
    }

  strcpy( tempModule->moduleName, moduleName );

  return tempModule;
  }


MODULEITEM *newModuleItem( MODULE *module )
  {
  MODULEITEM *tempModuleItem;

  tempModuleItem = malloc( sizeof( MODULEITEM ) );
  if( tempModuleItem == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newModuleItem()\n" );
    return NULL;
    }

  tempModuleItem->module = module;
  tempModuleItem->nextModuleItem = NULL;
 
  return tempModuleItem;
  }


MODULELIST *newModuleList( void )
  {
  MODULELIST *tempModuleList;

  tempModuleList = malloc( sizeof( MODULELIST ) );
  if( tempModuleList == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newModuleList()\n" );
    return NULL;
    }

  tempModuleList->firstModuleItem = NULL;

  return tempModuleList;
  }



/* -- String Methods ----------------------------------------------------------
*/

STRINGITEM *newStringItem( const char *string )
  {
  STRINGITEM *tempStringItem;

  tempStringItem = malloc( sizeof( STRINGITEM ) );
  if( tempStringItem == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newStringItem() [1]\n" );
    return NULL;
    }

  tempStringItem->string = malloc( strlen( string ) + 1 );
  if( tempStringItem == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newStringItem() [2]\n" );
    free( tempStringItem );
    return NULL;
    }

  strcpy( tempStringItem->string, string );
  tempStringItem->nextStringItem = NULL;

  return tempStringItem;
  }


void disposeStringItem( STRINGITEM *stringItem )
  {
  if( stringItem->string != NULL )
    free( stringItem->string );
  free( stringItem );
  }


void disposeStringItems( STRINGITEM *stringItem )
  {
  if( stringItem != NULL )
    {
    disposeStringItems( stringItem->nextStringItem );
    disposeStringItem( stringItem );
    }
  }


STRINGLIST *newStringList( void )
  {
  STRINGLIST *tempStringList;

  tempStringList = malloc( sizeof( STRINGLIST ) );
  if( tempStringList == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newStringList()\n" );
    return NULL;
    }

  tempStringList->firstStringItem = NULL;
  tempStringList->lastStringItem = NULL;

  return tempStringList;
  }


void disposeStringList( STRINGLIST *stringList )
  {
  disposeStringItems( stringList->firstStringItem );
  free( stringList );
  }


STRINGITEM *lookupStringItem( STRINGLIST *stringList, char *name )
  {
  STRINGITEM *item;

  item = stringList->firstStringItem;
  while( item != NULL )
    {
    if( strcmp( name, item->string ) == 0 )
      return item;

    item = item->nextStringItem;
    }

  /* nothing matched */
  return NULL;
  }

void addStringToList( const char *string, STRINGLIST *stringList )
  {
  STRINGITEM *stringItem;

  stringItem = newStringItem( string );
  if( stringItem == NULL )
    {
    fprintf( stderr, "newStringItem() returned NULL in addStringToList()\n" );
    return;
    }

  if( stringList->lastStringItem != NULL )
    {
    stringList->lastStringItem->nextStringItem = stringItem;
    stringList->lastStringItem = stringItem;
    }
  else
    {
    stringList->firstStringItem = stringItem;
    stringList->lastStringItem = stringItem;
    }
  }

char *newStringFromStringList( STRINGLIST *stringList )
  {
  size_t      totalStringLength;
  char       *newString;
  STRINGITEM *stringItem;

  totalStringLength = 0;

  stringItem = stringList->firstStringItem;
  while( stringItem != NULL )
    {
    totalStringLength += strlen( stringItem->string );
    stringItem = stringItem->nextStringItem;
    }
  
  newString = malloc( totalStringLength + 1 );
  if( newString == NULL )
    {
    fprintf( stderr, "Memory allocation failure in "
      "newStringFromStringList()\n" );
    return NULL;
    }
  
  newString[0] = '\0';  /* terminate string */

  stringItem = stringList->firstStringItem;
  while( stringItem != NULL )
    {
    strcat( newString, stringItem->string );
    stringItem = stringItem->nextStringItem;
    }

  return newString;
  }


BOOLEAN excluded( char *string, char *stringList[] )
  {
  char *listItem;
  int   i;

  i = 0;
  listItem = stringList[i];
  while( listItem != NULL )
    {
    if( strcmp( string, listItem ) == 0 )
      return True;
   
    listItem = stringList[i++];
    }

  return False;
  }


/* -- Module Methods ----------------------------------------------------------
*/

ROUTINE *lookupRoutine( ROUTINELIST *routineList, char *name )
  {
  ROUTINEITEM *item;

  item = routineList->firstRoutineItem;
  while( item != NULL )
    {
    if( strcmp( name, item->routine->routineName ) == 0 )
      return item->routine;

    item = item->nextRoutineItem;
    }

  /* nothing matched */
  return NULL;
  }


void addRoutineToList( ROUTINE *routine, ROUTINELIST *routineList )
  {
  ROUTINEITEM *routineItem;

  routineItem = newRoutineItem( routine );
  if( routineItem == NULL )
    {
    fprintf( stderr, "newRoutineItem() returned NULL in "
      "addRoutineToList()\n" );
    return;
    }

  routineItem->nextRoutineItem = routineList->firstRoutineItem;
  routineList->firstRoutineItem = routineItem;
  }


void addCallToList( CALLLIST *callList, char *name )
  {
  CALLITEM *callItem;
  CALL     *call;
#if 0
  int       size = 0;

  if( callList->truncated )
    return; /* List is already too large */
#endif

  callItem = callList->firstCallItem;
  while( callItem != NULL )
    {
/*        size++; */
      call = callItem->call;

      if( strcmp( name, call->name ) == 0 )
	return;  /* We already have this entry */

      callItem = callItem->nextCallItem;
    }

  /* If we got here then the call isn't in the list */

#if 0
  if( maxSize && size >= maxSize )
    {
      callList->truncated = True;
      return;
    }
#endif

  call = newCall( name );
  if( call == NULL )
    {
    fprintf( stderr, "Memory allocation failure in addCallToList() [1]\n" );
    return;
    }

  callItem = newCallItem( call );
  if( callItem == NULL )
    {
    fprintf( stderr, "Memory allocation failure in addCallToList() [2]\n" );
    disposeCall( call );
    return;
    }

  callItem->nextCallItem = callList->firstCallItem;
  callList->firstCallItem = callItem;
  }


void addCall( ROUTINELIST *routineList, ROUTINE *callerRoutine, char *callee )
  {
  ROUTINE *calleeRoutine;

  if( excluded( callee, callsExclusionList ) == False )
    addCallToList( callerRoutine->calls, callee );

  calleeRoutine = lookupRoutine( routineList, callee );
  if( calleeRoutine == NULL )
    {
    /* If there isn't already a routine here, then assume that this is
       a sort of forward reference, and create one. */

    calleeRoutine = newRoutine( callee );
    if( callee == NULL )
      {
      fprintf( stderr, "newRoutine() returned NULL in addCall() [1]\n" );
      return;
      }

    addRoutineToList( calleeRoutine, routineList );

    calleeRoutine = lookupRoutine( routineList, callee );
    if( calleeRoutine == NULL )
      {
      fprintf( stderr, "lookupRoutine() returned NULL in "
        "addCall() [2]\n" );
      }

    }
  addCallToList( calleeRoutine->calledFrom, callerRoutine->routineName );
  }


static FORTRANFILE *newFortranFile( char *fileName )
  {
  FORTRANFILE *tempFortranFile;

  tempFortranFile = malloc( sizeof( FORTRANFILE ) );
  if( tempFortranFile == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newFortranFile()\n" );
    return NULL;
    }

  tempFortranFile->file = newFileState( fileName );
  if( tempFortranFile->file == NULL )
    {
    fprintf( stderr, "newFileState() returned NULL in newFortranFile()\n" );
    free( tempFortranFile );
    return NULL;
    }

  tempFortranFile->currentLine = NULL;
  tempFortranFile->readPending = True;

  tempFortranFile->currentLineNumber = 0;
  tempFortranFile->currentByteOffset = 0;

  return tempFortranFile;
  }

void disposeFortranFile( FORTRANFILE *theFortranFile )
{
  /* discard currentLine? */
  disposeFileState( theFortranFile->file );
  free( theFortranFile );
}

void discardFortranLine( FORTRANFILE *file )
  {
  file->readPending = True;
  }


char *getFortranLine( FORTRANFILE *file )
  {
  char       *tempLine;
  char       *lineStart;
  char       *line;
  size_t      lineLength;
  char       *s;
  STRINGLIST *stringList;

  int         quoteState;
  const int   InsideQuote = 1;
  const int   OutsideQuote = 2;

  int         continuationState;
  const int   NotContinuation = 1;
  const int   Continuation = 2;

  if( file->readPending == True )
    {
    stringList = newStringList();
    if( stringList == NULL )
      {
      fprintf( stderr, "newStringList() returned NULL in getFortranLine()\n" );
      return NULL;
      }

    continuationState = NotContinuation;  /* Set initial state */
    line = getLine( file->file );
    if( line == NULL )
      return NULL;    /* Does this give a bogus line at the end of the file? */
    while( line != NULL )
      {
      lineLength = strlen( line );
      tempLine = malloc( lineLength + 1 );
      if( tempLine == NULL )
        {
        fprintf( stderr, "Memory allocation failure in "
          "getFortranLine()\n" );
        return NULL;
        }
      strcpy( tempLine, line );

      /* Get information about byte/line position */
/*
      file->currentLineNumber++;
      file->currentByteOffset += lineLength;
*/
      /* ... */

      /* Is the line a "C" comment? if so skip it */
      if( lineLength >= 7 && tempLine[0] == ' ' )
        {
        /* Is the line the first line or a continuation? */
        if( line[5] == ' ' && continuationState == Continuation )
          break;

        /* Deal with "!" comments */
        lineStart = tempLine + 6;
        s = lineStart;
        quoteState = OutsideQuote;
        while( *s != '\0' )
          { 
          if( *s == '!' )
            if( quoteState == OutsideQuote )
              {
              *s = '\0';  /* terminate line here */
              break;
              }

          if( *s == '\n' || *s == '\r' )
            {
            *s = '\0';  /* terminate line here */
            break;
            }

          if( *s == '\'' )
            {
            if( quoteState == OutsideQuote ) 
              quoteState = InsideQuote;
            else
              quoteState = OutsideQuote;
            }

          s++;
          }

        continuationState = Continuation;

        /* Since were still going add the results to the ... */
        addStringToList( lineStart, stringList );
        }

      free( tempLine );
 
      discardLine( file->file );
      file->currentLineNumber++;
      line = getLine( file->file );
      }

    if( file->currentLine != NULL )
      free( file->currentLine );

    file->currentLine = newStringFromStringList( stringList );
    file->readPending = False;

    disposeStringList( stringList );
    }

  return file->currentLine;
  }




/* -- Parsing Code ------------------------------------------------------------
*/

/*
 *  This also uppercase the string... 
 */
void mungeFortranLine( char *string )
  {
  int        quoteState;
  const int  InsideQuote  = 1;
  const int  OutsideQuote = 2;

  char      *source;
  char      *destination;

  destination = source = string;
  quoteState = OutsideQuote;
  while( *source != '\0' )
    {
    *destination = toupper( *source ); /* Uppercase in place */

    if( *source == '\'' )
      {
      if( quoteState == OutsideQuote )
        quoteState = InsideQuote;
      else
        quoteState = OutsideQuote;
      }
    else
      {
      if( quoteState == OutsideQuote )
        {
        *source = *destination;
        destination++;
        }
      }

    source++;
    }
  *destination = '\0';
  }


void processSubroutineCalls( char *line, ROUTINE *currentRoutine,
   ROUTINELIST *subroutineList )
  {
  PARSESTATE *parse;
  char       *name;

  parse = newParseState( line );
  if( parse == NULL )
    {
    fprintf( stderr, "newParseState() returned NULL in "
      "processSubroutineCalls()\n" );
    return;
    }

  while( empty( parse ) == False )
    {
    skipUntil( parse, "C" );
    if( prefix( parse, "CALL " ) == True )  /* What about tabs?? */
      {
      skipWord( parse, "CALL" );
      skip( parse, " \t\n\r" );
      name = extract( parse, "( " );
/*       addCall( routineList, currentRoutine->routineName, name ); */
      addCall( subroutineList, currentRoutine, name );
      free( name );
      }
    shiftChars( parse, 1 );
    }

  disposeParseState( parse );
  }

void processFunctionCalls( char *line, ROUTINE *currentRoutine,
   ROUTINELIST *functionList, STRINGLIST *stringList )
  {
  PARSESTATE *parse;
  PARSESTATE *subParse;
  char       *name;
  int         parenDepth;
  char        currentChar;

  parse = newParseState( line );
  if( parse == NULL )
    {
    fprintf( stderr, "newParseState() returned NULL in "
      "processFunctionCalls()\n" );
    return;
    }

  skip( parse, " \t\n\r" );  /* Whitespace */
  if( prefix( parse, "INTEGER" ) == True
     || prefix( parse, "REAL" ) == True
     || prefix( parse, "CHARACTER" ) == True
     || prefix( parse, "LOGICAL" ) == True
     || prefix( parse, "COMPLEX" ) == True )
    {
    /* make list of scalar variables and functions */
    skipUntil( parse, " \t\n\r" ); /* type */

    while( empty( parse ) == False )
      { 
      skip( parse, " \t\n\r" );  /* Whitespace */
      name = extract( parse, "*(, " );
      if( prefix( parse, "(" ) == False ) /* not an array */
	{
	if( lookupStringItem( stringList, name ) == NULL )
	  {
	  addStringToList( name, stringList );
	  }
	}
      free( name );
      skipUntil( parse, "," );
      shiftChars( parse, 1 ); /* ready for next scalar/function */
      }
    }
  else
    {
    skipUntil( parse, " =(" ); /* skip over keyword or = */
    shiftChars( parse, 1 );
    skip( parse, " (\t\n\r" );
    while( empty( parse ) == False )
      {
      name = extract( parse, " ()*/+-.,:\t\n\r" );
      if( prefix( parse, "(" ) == True )
	/* array, function call, or character substring */
	{
	/* check that it is not a character substring */
	shiftChars( parse, 1 ); /* over ( */
	subParse = newParseState( parse->currentPtr );
	parenDepth = 1;
	while( parenDepth > 0 )
	  {
	  skipUntil( subParse, "():," );
	  currentChar = *(subParse->currentPtr);
	  if( currentChar == '(' )
	    {
	    parenDepth++;	      
	    }
	  else if( parenDepth == 1 )
	    {
	    if( currentChar != ':' ) /* not a substring */
	      {
	      parenDepth = 0;	      
	      }
	    break;
	    }
	  else if( currentChar == ')' )
	    {
	    parenDepth--;	      
	    }
	  shiftChars( subParse, 1 );
	  }
	disposeParseState( subParse );
	if( parenDepth == 0 )  /* array or function call */
	  {
	  if( lookupStringItem( stringList, name ) != NULL )
	    {
	    addCall( functionList, currentRoutine, name );
	    }
	  }
        }
      free( name );
      skip( parse, " =()*/+-.,:\t\n\r" );
      }
    }
  disposeParseState( parse );
  }


void processExternalSymbols( char *line, ROUTINE *currentRoutine,
   ROUTINELIST *routineList )
  {
  PARSESTATE *parse;
  char       *name;

  parse = newParseState( line );
  skip( parse, " \t\n\r" );  /* Whitespace */
  if( parse == NULL )
    {
    fprintf( stderr, "newParseState() returned NULL in "
      "processExternalSymbols()\n" );
    return;
    }

  /* KAT 8Dec99: looking for externals also */
  if( prefix( parse, "EXTERNAL" ) == True )
    {
      skipWord( parse, "EXTERNAL" );
      skip( parse, " \t\n\r" );
      while( empty( parse ) == False )
	{
	  name = extract( parse, ", \t\n\r" );
	  addCall( routineList, currentRoutine, name );
	  skip( parse, ", \t\n\r" );
	}
    }
  disposeParseState( parse );
  }


void addRoutineDetailsToList( ROUTINELIST *routineList, char *routineName,
  char *moduleName )
  {
  ROUTINE *routine;

  routine = lookupRoutine( routineList, routineName );
  if( routine == NULL )
    {
    routine = newRoutine( routineName );
    if( routine == NULL )
      {
      fprintf( stderr, "newRoutine() returned NULL in "
        "addRoutineDetailsToList()\n" );
      return;
      }

    addRoutineToList( routine, routineList );
    }

  /* add details here */
  routine->moduleName = malloc( strlen( moduleName ) + 1 );
  if( routine->moduleName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in "
      "addRoutineDetailsToList()\n" );
    return;
    }

  strcpy( routine->moduleName, moduleName );
  }


void collectExternalSymbols( ROUTINELIST *externalList,
			     ROUTINELIST *subroutineList,
			     ROUTINELIST *functionList )
/* 8Dec99 KAT: collects called-from components in externalList and
   assembles them into subroutineList and functionList. */
{
  ROUTINEITEM *externalItem;
  ROUTINE     *routine;
  CALLITEM    *callItem;
  char        *externalName;

  externalItem = externalList->firstRoutineItem;
  while( externalItem != NULL )
    {
      externalName = externalItem->routine->routineName;
      if( ((routine = lookupRoutine( subroutineList, externalName ))
	   || (routine = lookupRoutine( functionList, externalName )))
	  == NULL )
	{
	  fprintf( stderr,
		   "CallTree Warning: No routine found for external %s.\n",
		   externalName );
	}
      else
	{
	  callItem = externalItem->routine->calledFrom->firstCallItem;
	  while( callItem != NULL )
	    {
	      addCallToList( routine->calledFrom, callItem->call->name );

	      callItem = callItem->nextCallItem;
	    }
	}
      externalItem = externalItem->nextRoutineItem;
    }
}


void dumpList( ROUTINELIST *routineList )
  {
  ROUTINEITEM *routineItem;
  ROUTINE     *routine;
  CALLITEM    *callItem;
  CALL        *call;

  routineItem = routineList->firstRoutineItem;
  while( routineItem != NULL )
    {
    routine = routineItem->routine;

    printf( "routine %s\n", routine->routineName );

    callItem = routine->calls->firstCallItem;
    if( callItem != NULL )
      {
      printf( "Calls:\n" );
      while( callItem != NULL )
        {
        call = callItem->call;

        printf( "  %s\n", call->name );

        callItem = callItem->nextCallItem;
        }
      }

    callItem = routine->calledFrom->firstCallItem;
    if( callItem != NULL )
      {
      printf( "Called From:\n" );
      while( callItem != NULL )
        {
        call = callItem->call;

        printf( "  %s\n", call->name );

        callItem = callItem->nextCallItem;
        }
      }

    printf( "\n" );

    routineItem = routineItem->nextRoutineItem;
    }
  }


void parseFile( ROUTINELIST *subroutineList, ROUTINELIST *functionList,
		ROUTINELIST *externalList, char *fileName )
  {
  FORTRANFILE *file;
  char        *line;
  PARSESTATE  *parse;
  char        *routineName;
  ROUTINE     *currentRoutine;
  char        *moduleName;
  STRINGLIST  *stringList;

  int fortranState;
  /* Damn C 
  const int OutsideRoutine = 1;
  const int InsideSubroutine = 2;
  const int InsideFunction = 3;
  */
  
  moduleName = malloc( strlen( fileName ) + 1 );  /* Overallocates */
  if( moduleName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in parseFile()\n" );
    return;
    }
  else
    {
    char *dir, *slash, *start, *end;
    char *source, *destination;
    const char src_dir[] = "source";
    const int src_dir_len = strlen( src_dir );

    /*
      To determine a module name use the subpath of the file relative to
      the last directory that looks like the source directory.
    */

    start = NULL;
    for( dir = fileName;
	 slash = strchr( dir, '/' ); /* UNIX specific pathname */
	 )
      {
	if( slash - dir == src_dir_len
	    && 0 == strncmp( dir, src_dir, src_dir_len ) )
	  {
	    start = slash + 1;
	  }

	dir = slash + 1;
      }

    if( start == NULL )
      {
	if( *start == '/' )
	  { /* Absolute path: just use the tail of the filename */
	    fprintf( stderr,
		     "Can't determine module from filename %s in parseFile()\n",
		     fileName );
	    start = strrchr( start, '/') + 1;
	  }
	else
	  { /* relative path; use the whole subpath */
	    start = fileName;
	  }
      }
    
    /* Trim the tail unless there is no directory component */
    end = strrchr( start, '/' );
    if( end == NULL ) end = start + strlen(start);

    for( source = start, destination = moduleName;
	 source < end;
	 source++, destination++ )
      {
	*destination = toupper( *source );
      }

    *destination = '\0';
    }

  file = newFortranFile( fileName );
  if( file == NULL )
    {
    fprintf( stderr, "newFortranFile() returned NULL in parseFile()\n" );
    /* Could return error if this routine had a return code */
    return;
    }

  fortranState = OutsideRoutine;
  line = getFortranLine( file );
  while( line != NULL )
    {
    mungeFortranLine( line );

    parse = newParseState( line );
    if( parse == NULL )
      {
      fprintf( stderr, "newParseState() returned NULL in "
        "parseFile()\n" );
      return;
      }

    switch( fortranState )
      {
      case OutsideRoutine:
        skip( parse, " \t\n\r" );  /* Whitespace */
        if( prefix( parse, "SUBROUTINE" ) == True )
          {
          skipWord( parse, "SUBROUTINE" );
          skip( parse, " \t\n\r" );  /* Whitespace */
          routineName = extract( parse, "( " );
          addRoutineDetailsToList( subroutineList, routineName, moduleName );
          currentRoutine = lookupRoutine( subroutineList, routineName );
          if( currentRoutine == NULL )
            {
            fprintf( stderr, "lookupRoutine() returned NULL in "
              "parseFile()\n" );
            return;
            }
          free( routineName ); /* Hmmm */
	  stringList = newStringList();
	  if( stringList == NULL )
	    {
	    fprintf( stderr, "newStringList() returned NULL in parseFile()\n" );
	    return;
	    }
          fortranState = InsideSubroutine;
          }

	skipUntil( parse, " \t\n\r" ); /* Type */
        skip( parse, " \t\n\r" );  /* Whitespace */
        if( prefix( parse, "FUNCTION" ) == True )
          {
          skipWord( parse, "FUNCTION" );
          skip( parse, " \t\n\r" );  /* Whitespace */
          routineName = extract( parse, "( " );
          addRoutineDetailsToList( functionList, routineName, moduleName );
          currentRoutine = lookupRoutine( functionList, routineName );
          if( currentRoutine == NULL )
            {
            fprintf( stderr, "lookupRoutine() returned NULL in "
              "parseFile()\n" );
            return;
            }
          free( routineName ); /* Hmmm */
	  stringList = newStringList();
	  if( stringList == NULL )
	    {
	    fprintf( stderr, "newStringList() returned NULL in parseFile()\n" );
	    return;
	    }
          fortranState = InsideFunction;
          }

        break;

      case InsideFunction:
      case InsideSubroutine:
        skip( parse, " \t\n\r" );
        if( prefix( parse, "END" ) == True )
          {
          skipWord( parse, "END" );
          skip( parse, " \t\n\r" );
          if( empty( parse ) == True )
	    {
	    disposeStringList( stringList );
            fortranState = OutsideRoutine;
	    }
          }
	else
	  {
	  processSubroutineCalls( line, currentRoutine, subroutineList );
	  processFunctionCalls( line, currentRoutine, functionList,
			       stringList );
	  processExternalSymbols( line, currentRoutine, externalList );
	  }
        break;

      default:
        fprintf( stderr, "Illegal program state in parseFile()\n" );
        return;
      }

    disposeParseState( parse );
      
    discardFortranLine( file );
    line = getFortranLine( file );
    }

  disposeFortranFile( file );
  if( fortranState != OutsideRoutine )
    {
    disposeStringList( stringList );
    }

  free( moduleName );   /* loose free - messy */
  }



void parseSourceFiles( ROUTINELIST *subroutineList, ROUTINELIST *functionList,
		      char *moduleListFileName )
  {
  FILESTATE   *moduleList;
  ROUTINELIST *externalList;
  PARSESTATE  *parse;
  char        *line;
  char        *fileName;

  moduleList = newFileState( moduleListFileName );
  if( moduleList == NULL )
    {
    fprintf( stderr, "newFileState() returned NULL in parseSourceFiles()\n" );
    return;
    }

  externalList = newRoutineList();

  line = getLine( moduleList );
  while( line != NULL )
    {
    parse = newParseState( line );
    if( parse == NULL )
      {
      fprintf( stderr, "newParseState() returned NULL in "
        "parseSourceFiles()\n" );
      return;
      }

    skip( parse, " \t\n\r" );
    fileName = extract( parse, " \t\n\r" );
    parseFile( subroutineList, functionList, externalList, fileName );
    disposeParseState( parse );

    discardLine( moduleList );
    line = getLine( moduleList );
    }

  disposeFileState( moduleList );

  collectExternalSymbols( externalList, subroutineList, functionList );

  disposeRoutineList( externalList );
  }



/*
 *  This routine gives the term "inefficient" new meaning, in a program
 *  that is appallingly slow as is...
 */
TAG *newCallListTag( const char *tagName, CALLLIST *callList )
  {
  TAG      *tempTag;
  CALLITEM *callItem;
  CALL     *call;

  int       callsListed;
  int       stringLength;
  char     *listString;
  int       lineCalls;
  int       currentStringLength;
  int       currentLineLength;

  char     *source;
  char     *destination;

  const int MaximumWidth = 60;
    /* --  maximum length of line */
/*    char truncation[] = ", ..."; */

  if( callList->firstCallItem == NULL )
    return NULL;         /* No tag is generated */

  /* find length of output string, and number of output calls */
  callItem = callList->firstCallItem;
  callsListed = 0;
  stringLength = 0;
  lineCalls = 0;
  currentLineLength = 0;
  while( callItem != NULL )
    {
    call = callItem->call;

    callsListed++;

    currentStringLength = strlen( call->name );

    if( callsListed > 1 )
      currentStringLength += 2;  /* <comma space> needed */
      
    if( lineCalls > 1 )
      if( currentLineLength + currentStringLength > MaximumWidth )
	{
          stringLength += 1; /* insert a newline */
          currentLineLength = 0;
          lineCalls = 0;
	}

    currentLineLength += currentStringLength;
    stringLength += currentStringLength;

    lineCalls++;

    callItem = callItem->nextCallItem;
    }

#if 0
  if( callList->truncated == True )
    {
      stringLength += strlen( truncation );
    }
#endif

  if( callsListed == 0 )
    return NULL;         /* No tag is generated */

  /* generate string */
  listString = malloc( stringLength + 1 );
  if( listString == NULL )
    {
    fprintf( stderr, "Memory allcoation failure in newCallListTag()\n" );
    return NULL;
    }

  callItem = callList->firstCallItem;
  callsListed = 0;
  lineCalls = 0;
  currentLineLength = 0;
  destination = listString;
  while( callItem != NULL )
    {
    call = callItem->call;

    callsListed++;

    currentStringLength = strlen( call->name );

    if( callsListed > 1 )
      {
        *destination++ = ',';  /* comma */
        *destination++ = ' ';  /* space */
        currentStringLength += 2;
      }

    if( lineCalls > 1 )
      if( currentLineLength + currentStringLength > MaximumWidth )
	{
          *destination++ = '\n'; /* insert a newline */
          currentLineLength = 0;
          lineCalls = 0;
	}

    source = call->name;
    while( *source )
      *destination++ = *source++;

    currentLineLength += currentStringLength;

    lineCalls++;

    callItem = callItem->nextCallItem;
    }
#if 0
  if( callList->truncated == True )
    {
      source = truncation;
      while( *source )
	*destination++ = *source++;      
    }
#endif
  *destination = '\0';

  /* Assertion */
  if( stringLength != strlen( listString ) )
    fprintf( stderr, "Assertion failed in newCallListTag() - "
      "stringLength != strlen( listString )\n" );

  tempTag = newTag();
  if( tempTag == NULL )
    {
    fprintf( stderr, "newTag() returned NULL in newCallListTag()\n" );
    free( listString );
    return NULL;
    }
  
  tempTag->tag = malloc( strlen( tagName ) + 1 );
  if( tempTag->tag == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newCallListTag() [2]\n" );
    free( listString );
    disposeTag( tempTag );
    return NULL;
    }

  strcpy( tempTag->tag, tagName );
  tempTag->body = listString;

  return tempTag;
  }


int lexCompareCallItems( CALLITEM **firstItem, CALLITEM **secondItem )
  {
  return strcmp( (*firstItem)->call->name, (*secondItem)->call->name );
  }


void sortCallList( CALLLIST *callList )
  {
  CALLITEM  *callItem;
  CALLITEM **callItemArray;
  int        numCallItems;
  int        i;

  /* step 1. count the number of elements */
  numCallItems = 0;
  callItem = callList->firstCallItem;
  while( callItem != NULL )
    {
    numCallItems++;
    callItem = callItem->nextCallItem;
    }

  if( numCallItems < 2 )
    return;  /* don't sort zero or one items */

  /* step 2. create array for sorting */
  callItemArray = malloc( sizeof( CALLITEM * ) * numCallItems );
  if( callItemArray == NULL )
    {
    fprintf( stderr, "Memory allocation failure in sortCallList()\n" );
    return;
    }

  /* step 3. dump list contents to array */
  i = 0;
  callItem = callList->firstCallItem;
  while( callItem != NULL )
    {
    callItemArray[ i++ ] = callItem;
    callItem = callItem->nextCallItem;
    }

  /* step 4. invoke qsort on array */
  qsort( callItemArray, numCallItems, sizeof( CALLITEM * ), 
    (int (*)(const void *, const void *)) lexCompareCallItems );

  /* step 5. re-thread list from array */
  callList->firstCallItem = callItemArray[0];
  for( i = 0; i < numCallItems - 1; i++ )
    callItemArray[i]->nextCallItem = callItemArray[i+1];
  callItemArray[numCallItems - 1]->nextCallItem = NULL;

  /* step 6. free array */
  free( callItemArray );
  }


TAG *newDefinedInTag( char *moduleName )
  {
  TAG *tempTag;

  const char *DefinedInTagName = "Defined-in";
  
  tempTag = newTagFromContents( DefinedInTagName, moduleName );
  if( tempTag == NULL )
    {
    fprintf( stderr,
	     "newTagFromContents() returned NULL in newDefinedInTag()\n" );
    return NULL;
    }

  return tempTag;
  }


void outputCallTreeInformation( ROUTINELIST *routineList, 
  char *routineFileName, REDUCTION_FUNCTION *reduce, char *groupTag )
  {
#ifndef WORKING
  /* Want to open a collection, and add information to it... */
  COLLECTION  *collection;
  GROUPLIST   *groupList;
  GROUP       *group;
  ROUTINEITEM *routineItem;
  ROUTINE     *routine;
  TAG         *callsTag;
  TAG         *calledFromTag;
  TAG         *definedInTag;

  collection = newCollectionFromFile( routineFileName );
  if( collection == NULL )
    {
    fprintf( stderr, "newCollectionFromFile() returned NULL in "
      "outputCallTreeInformation()\n" );
    return;
    }

  groupList = lookupGroupList( collection, groupTag );
  if( groupList == NULL )
    {
    fprintf( stderr, "lookupGroupList() returned NULL in "
      "outputCallTreeInformation()\n" );
    disposeCollection( collection );
    return;
    }

  /* Go through each routine and add a tag to the corresponding
     variable group for "Calls:" and "Called-From:" (if they exist) */
  routineItem = routineList->firstRoutineItem;
  while( routineItem != NULL )
    {
    routine = routineItem->routine;

    group = lookupReducedGroup( groupList, routine->routineName, reduce );

    if( group == NULL )
      {
	TAG *tag;

	tag = newTagFromContents( groupTag, routine->routineName );
	if( tag == NULL )
	  {
	    fprintf( stderr,
 "newTagFromContents() returned NULL in outputCallTreeInformation()\n" );
	    disposeCollection( collection );
	    return;
	  }

	group = newGroupFromHeadTag( tag, groupTag );
	if( group == NULL )
	  {
	    fprintf( stderr,
 "newGroupFromHeadTag() returned NULL in outputCallTreeInformation()" );
	    disposeCollection( collection );
	    disposeTag( tag );
	    return;
	  }

	addGroupToGroupList( group, groupList );

	if( routine->moduleName != NULL )
	  fprintf( stderr,
	    "CallTree Warning: No documentation found for routine %s.\n",
	    routine->routineName );
      }
#if 0
      fprintf( stderr,
        "CallTree Warning: trying to add call data to the routine %s\n"
        "but no such routine exists in the database. Please check your\n"
        "source files for appropriate documentation lines.\n",
        routine->routineName );
      }
    else
      {
#endif
      sortCallList( routine->calls );
      callsTag = newCallListTag( "Calls", routine->calls );
      sortCallList( routine->calledFrom );
      calledFromTag = newCallListTag( "Called-From", routine->calledFrom );

      if( callsTag != NULL )
        addTagToGroup( callsTag, group );
                                 /* multiple return semantics - NULL means 
                                    either a failure (as in a malloc failed)
                                    *OR* that there was no data to return.
                                    Basically untidy - fix later */
    
      if( calledFromTag != NULL )
        addTagToGroup( calledFromTag, group );  /* Ditto */

      if( routine->moduleName != NULL )
        {
        definedInTag = newDefinedInTag( routine->moduleName );
        if( definedInTag != NULL )
          addTagToGroup( definedInTag, group ); 
        }
/*        } */

    routineItem = routineItem->nextRoutineItem;
    }

  writeCollection( collection, routineFileName );
#else
  /* Dumps list to stdout */
  dumpList( routineList );
#endif
  }


MODULE *lookupModule( char *moduleName, MODULELIST *moduleList )
  {
  MODULEITEM *moduleItem;

  moduleItem = moduleList->firstModuleItem;
  while( moduleItem != NULL )
    {
    if( strcmp( moduleName, moduleItem->module->moduleName ) == StringMatch )
      {
      return moduleItem->module;
      }

    moduleItem = moduleItem->nextModuleItem;
    }

  /* No match found */
  return NULL;
  }


void addModuleToModuleList( MODULE *module, MODULELIST *moduleList )
  {
  MODULEITEM *moduleItem;

  moduleItem = newModuleItem( module );
  if( moduleItem == NULL )
    {
    fprintf( stderr, "newModuleItem() returned NULL in "
      "addModuleToModuleList()\n" );
    return;
    }

  moduleItem->nextModuleItem = moduleList->firstModuleItem;
  moduleList->firstModuleItem = moduleItem;
  }


void addRoutineToModule( ROUTINE *routine, MODULE *module )
  {
  addRoutineToList( routine, module->routineList );
  }


void addRoutineToModuleList( ROUTINE *routine, MODULELIST *moduleList )
  {
  MODULE *module;

  module = lookupModule( routine->moduleName, moduleList );
  if( module == NULL )
    {
    module = newModule( routine->moduleName );
    if( module == NULL )
      {
      fprintf( stderr, "newModule() returned NULL in "
        "addRoutineToModuleList()\n" );
      return;
      }
    addModuleToModuleList( module, moduleList );
    }

  addRoutineToModule( routine, module );
  }


int lexCompareRoutineItems( ROUTINEITEM **firstItem, ROUTINEITEM **secondItem )
  {
  return strcmp( (*firstItem)->routine->routineName, 
    (*secondItem)->routine->routineName );
  }


void sortRoutineList( ROUTINELIST *routineList )
  {
  ROUTINEITEM  *routineItem;
  ROUTINEITEM **routineItemArray;
  int           numRoutineItems;
  int           i;

  /* step 1. count the number of elements */
  numRoutineItems = 0;
  routineItem = routineList->firstRoutineItem;
  while( routineItem != NULL )
    {
    numRoutineItems++;
    routineItem = routineItem->nextRoutineItem;
    }

  if( numRoutineItems < 2 )
    return;  /* don't sort zero or one items */

  /* step 2. create array for sorting */
  routineItemArray = malloc( sizeof( ROUTINEITEM * ) * numRoutineItems );
  if( routineItemArray == NULL )
    {
    fprintf( stderr, "Memory allocation failure in sortRoutineList()\n" );
    return;
    }

  /* step 3. dump list contents to array */
  i = 0;
  routineItem = routineList->firstRoutineItem;
  while( routineItem != NULL )
    {
    routineItemArray[ i++ ] = routineItem;
    routineItem = routineItem->nextRoutineItem;
    }

  /* step 4. invoke qsort on array */
  qsort( routineItemArray, numRoutineItems, sizeof( ROUTINEITEM * ), 
    (int (*)(const void *, const void *)) lexCompareRoutineItems );

  /* step 5. re-thread list from array */
  routineList->firstRoutineItem = routineItemArray[0];
  for( i = 0; i < numRoutineItems - 1; i++ )
    routineItemArray[i]->nextRoutineItem = routineItemArray[i+1];
  routineItemArray[numRoutineItems - 1]->nextRoutineItem = NULL;

  /* step 6. free array */
  free( routineItemArray );
  }


/* Great. Cargo cult programming. And I wasn't to pleased with the 
   original function to be honest... */
/* To fix this create a datatype that contains the output string length,
   the output string and the nessessary state, then call a method
   like "addNameToList( name, nameList );". When done call
   "getString( nameList );" */
TAG *newRoutinesTag( ROUTINELIST *routineList )
  {
  TAG         *tempTag;
  ROUTINEITEM *routineItem;
  ROUTINE     *routine;

  int       routinesListed;
  int       stringLength;
  char     *listString;
  int       lineRoutines;
  int       currentStringLength;
  int       currentLineLength;

  char     *source;
  char     *destination;

  const int MaximumWidth = 60;
    /* --  maximum length of line */

  const char *RoutinesTagName = "ROUTINES";

  if( routineList->firstRoutineItem == NULL )
    return NULL;         /* No tag is generated */

  /* find length of output string, and number of output routines */
  routineItem = routineList->firstRoutineItem;
  routinesListed = 0;
  stringLength = 0;
  lineRoutines = 0;
  currentLineLength = 0;
  while( routineItem != NULL )
    {
    routine = routineItem->routine;

    routinesListed++;

    currentStringLength = strlen( routine->routineName );

    if( routinesListed > 1 )
      currentStringLength += 2;  /* <comma space> needed */
    
    if( lineRoutines > 1 )
      if( currentLineLength + currentStringLength > MaximumWidth )
        {
        stringLength += 1; /* insert a newline */
        currentLineLength = 0;
        lineRoutines = 0;
        }

    currentLineLength += currentStringLength;
    stringLength += currentStringLength;

    lineRoutines++;

    routineItem = routineItem->nextRoutineItem;
    }

  if( routinesListed == 0 )
    return NULL;         /* No tag is generated */

  /* generate string */
  listString = malloc( stringLength + 1 );
  if( listString == NULL )
    {
    fprintf( stderr, "Memory allcoation failure in newRoutinesTag()\n" );
    return NULL;
    }

  routineItem = routineList->firstRoutineItem;
  routinesListed = 0;
  lineRoutines = 0;
  currentLineLength = 0;
  destination = listString;
  while( routineItem != NULL )
    {
    routine = routineItem->routine;

    routinesListed++;

    currentStringLength = strlen( routine->routineName );

    if( routinesListed > 1 )
      {
      *destination++ = ',';  /* comma */
      *destination++ = ' ';  /* space */
      currentStringLength += 2;
      }

    if( lineRoutines > 1 )
      if( currentLineLength + currentStringLength > MaximumWidth )
        {
        *destination++ = '\n'; /* insert a newline */
        currentLineLength = 0;
        lineRoutines = 0;
        }

    source = routine->routineName;
    while( *source )
      *destination++ = *source++;

    currentLineLength += currentStringLength;

    lineRoutines++;

    routineItem = routineItem->nextRoutineItem;
    }
  *destination = '\0';

  /* Assertion */
  if( stringLength != strlen( listString ) )
    fprintf( stderr, "Assertion failed in newRoutinesTag() - "
      "stringLength != strlen( listString )\n" );

  tempTag = newTag();
  if( tempTag == NULL )
    {
    fprintf( stderr, "newTag() returned NULL in newRoutinesTag()\n" );
    free( listString );
    return NULL;
    }
  
  tempTag->tag = malloc( strlen( RoutinesTagName ) + 1 );
  if( tempTag->tag == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newRoutinesTag() [2]\n" );
    free( listString );
    disposeTag( tempTag );
    return NULL;
    }

  strcpy( tempTag->tag, RoutinesTagName );
  tempTag->body = listString;

  return tempTag;
  }


#define FUNCTION "writeModuleInformationToCollection"
void writeModuleInformationToCollection( MODULELIST *moduleList, 
  char *moduleBaseName )
  {
  COLLECTION  *collection;
  GROUPLIST   *groupList;
  GROUP       *group;
  MODULEITEM  *moduleItem;
  MODULE      *module;
  TAG         *routinesTag;

  collection = newCollectionFromFile( moduleBaseName );
  if( collection == NULL )
    {
    fprintf( stderr, "newCollectionFromFile() returned NULL in "
      "writeModuleInformationToCollection()\n" );
    return;
    }

  groupList = lookupGroupList( collection, "MODULE" );
  if( groupList == NULL )
    {
    fprintf( stderr, "lookupGroupList() returned NULL in "
      "writeModuleInformationToCollection()\n" );
    disposeCollection( collection );
    return;
    }

  moduleItem = moduleList->firstModuleItem;
  while( moduleItem != NULL )
    {
    module = moduleItem->module;

    group = lookupReducedGroup( groupList, module->moduleName, reduceModule );
    if( group == NULL )
      {
	TAG *tag;

	tag = newTagFromContents( "MODULE", module->moduleName );
	if( tag == NULL )
	  {
	    fprintf( stderr,
 "newTagFromContents() returned NULL in " FUNCTION "\n" );
	    disposeCollection( collection );
	    return;
	  }

	group = newGroupFromHeadTag( tag, "MODULE" );
	if( group == NULL )
	  {
	    fprintf( stderr,
 "newGroupFromHeadTag() returned NULL in " FUNCTION "\n" );
	    disposeTag( tag );
	    disposeCollection( collection );
	    return;
	  }

	addGroupToGroupList( group, groupList );

	fprintf( stderr,
		 "CallTree Warning: No documentation found for module %s.\n",
		 module->moduleName );

      }
#if 0
      fprintf( stderr,
        "CallTree Warning: trying to add subroutine data to the module %s\n"
        "but no such module exists in the database. Please check your\n"
        "source files for appropriate documentation lines.\n",
        module->moduleName );
      }
    else
      {
#endif
      {
      sortRoutineList( module->routineList );
      routinesTag = newRoutinesTag( module->routineList );

      if( routinesTag != NULL )
        addTagToGroup( routinesTag, group );
      }

    moduleItem = moduleItem->nextModuleItem;
    }

  writeCollection( collection, moduleBaseName );
  }
#undef FUNCTION


void outputModuleInformation( ROUTINELIST *subroutineList,
			      ROUTINELIST *functionList,
			      char *moduleBaseName )
  {
  MODULELIST  *moduleList;
  ROUTINEITEM *routineItem;
  ROUTINELIST *routineLists[2];
  int i;

  moduleList = newModuleList();
  if( moduleList == NULL )
    {
    fprintf( stderr, "newModuleList() returned NULL in "
      "outputModuleInformation()\n" ); 
    return;
    }

  routineLists[0] = functionList;
  routineLists[1] = subroutineList;

  for( i = 0; i < 2; i++ )
    {
      routineItem = routineLists[i]->firstRoutineItem;

      while( routineItem != NULL )
	{
	  /* Only add routines with Module information */
	  if( routineItem->routine->moduleName != NULL )
	    addRoutineToModuleList( routineItem->routine, moduleList );

	  routineItem = routineItem->nextRoutineItem;
	}
    }

  writeModuleInformationToCollection( moduleList, moduleBaseName );
  }


void outputMasterModuleList( char *moduleBaseName, char *listFile )
{
  COLLECTION  *collection;
  GROUPLIST *groupList;
  GROUPITEM *groupItem; 
  GROUP     *group;
  TAG       *tag;
  int       counter,i,total;
  char      *moduleName, *reducedName;
  FILE      *file;
  const char *comment;

  file = fopen( listFile, "wt" );
  if( file == NULL )
  {
    fprintf( stderr, "Failed to open file %s in outputMasterModuleList(): %s\n",
	     listFile, strerror(errno) );
    return;
  }

  fprintf( file, "<HTML><HEAD><TITLE>CMISS Module List</TITLE></HEAD>" );
  fprintf( file, "<BODY><H1>CMISS Module List</H1><HR><P>" );
  fprintf( file, "<UL>" );

  collection = newCollectionFromFile( moduleBaseName );
  if( collection == NULL )
    {
    fprintf( stderr, "newCollectionFromFile() returned NULL in "
      "writeModuleInformationToCollection()\n" );
    return;
    }

  groupList = lookupGroupList( collection, "MODULE" );
  if( groupList == NULL )
  {
    fprintf( stderr, "lookupGroupList() failed in outputGroupListFile()\n" );
    return;
  }


  total=0;
  groupItem = groupList->firstGroupItem;
  while( groupItem != NULL )
  {
    total++;
    groupItem = groupItem->nextGroupItem;
  }
  
  groupItem = groupList->firstGroupItem;
  counter=total;
  i=1;
  while( counter != 0 )
  {
    while( counter != 1 )
    {
      counter--;
      groupItem = groupItem->nextGroupItem;
    }
    
    group = groupItem->group;

    tag = lookupTag( group, "MODULE" );
    moduleName = tag->body;

    reducedName = malloc( strlen( moduleName ) +1 );
    reduceModule( reducedName, moduleName);

    fprintf( file, "<LI><A HREF=\"%s?name=%s\">%s</A>\n",
      LOOKUP_PROGRAM_URLPATH, reducedName, moduleName );
    /* -- NB, It would be better if we urlencoded the names first... */

    tag = lookupTag( group, "DESCRIPTION" );
    if( tag != NULL )
      {
	comment = tag->body;

	fprintf( file, "- %s\n", comment );
      }

    counter=total-i;
    i++;
    groupItem = groupList->firstGroupItem;
  }

#if 0
  /* CS new Mon Jan  4 18:29:39 "NDT 1999 */
  fprintf( file, "<LI><A HREF=\"%s?name=cmiss_archive1.f\">CMISS_ARCHIVE1</A> - Contains archived code from modules FE00 -> FE09",LOOKUP_PROGRAM_URLPATH );
  
  fprintf( file, "<LI><A HREF=\"%s?name=cmiss_archive2.f\">CMISS_ARCHIVE2</A> - Contains archived code from modules FE10 -> FE19",LOOKUP_PROGRAM_URLPATH );

  fprintf( file, "<LI><A HREF=\"%s?name=cmiss_archive3.f\">CMISS_ARCHIVE3</A> - Contains archived code from modules FE20 -> FE29",LOOKUP_PROGRAM_URLPATH );
  
  fprintf( file, "<LI><A HREF=\"%s?name=cmiss_archive4.f\">CMISS_ARCHIVE4</A> - Contains archived code from modules FE30 -> FEUSER",LOOKUP_PROGRAM_URLPATH );
#endif  

  fprintf( file, "</UL>" );
}


void printHelp( char *programName )
  {
  printf( "Syntax: %s [-s <subroutine list base name>] "
    "[-m <module list base name>] [-l <module list filename>]\n",
    programName );
  }



/* -- Program Entry Point -----------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
  char *subroutineFileName = DefaultSubroutineFileName;
  char *functionFileName = DefaultFunctionFileName;
  char *moduleListFileName = DefaultModuleListFileName;
  char *moduleBaseName = DefaultModuleBaseName;

  ROUTINELIST *subroutineList;
  ROUTINELIST *functionList;

  int   arg;
  char *optPtr;

  if( argc > 1 )
    {
    arg = 1;

    while( arg < argc )
      {
      if( argv[ arg ] != NULL && argv[ arg ][0] == '-' )
        {
        /* an option */
        optPtr = &(argv[ arg ][1]);

        while( *optPtr != NULL )
          {
          switch( *optPtr )
            {
            case 's': 
              arg++;
              if( argv[ arg ] != NULL )
                {
                subroutineFileName = argv[ arg ];
                }
              else
                {
                fprintf( stderr, "Error: No Parameter for -s option\n" );
                printHelp( argv[0] );
                exit( 1 );
                }
              break;

            case 'f': 
              arg++;
              if( argv[ arg ] != NULL )
                {
                functionFileName = argv[ arg ];
                }
              else
                {
                fprintf( stderr, "Error: No Parameter for -f option\n" );
                printHelp( argv[0] );
                exit( 1 );
                }
              break;

            case 'm': 
              arg++;
              if( argv[ arg ] != NULL )
                {
                moduleBaseName = argv[ arg ];
                }
              else
                {
                fprintf( stderr, "Error: No Parameter for -m option\n" );
                printHelp( argv[0] );
                exit( 1 );
                }
              break;

            case 'l': 
              arg++;
              if( argv[ arg ] != NULL )
                {
                moduleListFileName = argv[ arg ];
                }
              else
                {
                fprintf( stderr, "Error: No Parameter for -l option\n" );
                printHelp( argv[0] );
                exit( 1 );
                }
              break;

            case 'h':
              printHelp( argv[0] );
              return 0;

            default: 
              fprintf( stderr, "Illegal Option: %c\n", *optPtr );
              printHelp( argv[0] );
              exit( 1 );
            }

          optPtr++;
          }
        }
      else
        {
        if( argv[ arg ] != NULL )
          {
          fprintf( stderr, "Unknown Parameter: %s\n", argv[ arg ] );
          printHelp( argv[0] );
          exit( 1 );
          }
        }

      arg++;
      }
    }

  subroutineList = newRoutineList();
  functionList = newRoutineList();
  if( subroutineList == NULL || functionList == NULL )
    {
    fprintf( stderr, "newRoutineList() returnd NULL in main()\n" );
    exit( 1 );
    }

  parseSourceFiles( subroutineList, functionList, moduleListFileName );
  outputCallTreeInformation( subroutineList, subroutineFileName,
			    reduceSubroutine, "SUBROUTINE" );
  outputCallTreeInformation( functionList, functionFileName,
			    reduceFunction, "FUNCTION" );
  outputModuleInformation( subroutineList, functionList, moduleBaseName );
  outputMasterModuleList( moduleBaseName, MasterModuleListFile );

  return 0;
  }
