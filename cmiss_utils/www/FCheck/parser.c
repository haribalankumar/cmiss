/*
 *  Parser.c
 *
 *    Code to parse the cmiss check and common block files
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf() */

#include <stdlib.h>
    /* For: malloc(), free(), NULL */

#include <string.h>
    /* For: strcpy(), strlen(), strcpy(), etc */

#include <time.h>
    /* For: time_t, struct tm, time(), localtime(), strftime() */

#include "lineio.h"
    /* For: FILESTATE, newFileState(), getLine(), discardLine() */

#include "simple-parser.h"
    /* For: PARSESTATE, newParseState(), prefix(), skip(), extract() etc */



/* -- Module Constants --------------------------------------------------------
*/

char *WhiteSpace = " \t\n\r";
char *SetNeverUsed = "  Elements set but never used:";
char *UsedNeverSet = "  Elements used but never set:";
char *NeverUsedNeverSet = "  Elements never used, never set:";



/* -- Private Module Datatypes ------------------------------------------------
*/ 

typedef struct STRINGITEM
  {
  char *string;

  struct STRINGITEM *nextString;
  }
STRINGITEM;


typedef struct
  {
  STRINGITEM *firstString;
  }
STRINGLIST;


typedef struct ROUTINE
  {
  char *name;
  char *sourceFile;

  struct ROUTINE *nextRoutine;
  }
ROUTINE;


typedef struct COMMON
  {
  char       *name;
  BOOLEAN     globallyUnused;
  ROUTINE    *firstUnusedRoutine;
  STRINGLIST *variablesSetNeverUsed;
  STRINGLIST *variablesUsedNeverSet;
  STRINGLIST *variablesNeverUsedNeverSet;

  struct COMMON *nextCommonBlock;
  }
COMMON;


typedef struct INCLUDE
  {
  char    *name;
  char    *fileName;
  COMMON  *firstCommonBlock;
  ROUTINE *firstUnusedRoutine;

  struct INCLUDE *nextInclude;
  }
INCLUDE;


typedef struct INCLUDENAME
  {
  char *name;
  
  struct INCLUDENAME *nextInclude;
  }
INCLUDENAME;


typedef struct ROUTINENAME
  {
  char        *name;
  char        *file;
  INCLUDENAME *firstInclude;

  struct ROUTINENAME *nextRoutine;
  }
ROUTINENAME;


typedef struct
  {
  INCLUDE     *firstInclude;
  ROUTINENAME *firstInvertedRoutine;
  }
INCLUDEITEMS;



/* -- Constructor, Destructors and Duplicators --------------------------------
*/

STRINGITEM *newStringItem( char *string )
  {
  STRINGITEM *tempStringItem;

  tempStringItem = malloc( sizeof( STRINGITEM ) );
  if( tempStringItem == NULL )
    {
    fprintf( stderr, "Memory allocation error in newStringItem()\n" );
    return NULL;
    }

  tempStringItem->string = string;
  tempStringItem->nextString = NULL;

  return tempStringItem;
  }


void disposeStringItem( STRINGITEM *theStringItem )
  {
  free( theStringItem->nextString );
  free( theStringItem );
  }


void disposeStringItemList( STRINGITEM *theItem )
  {
  if( theItem != NULL )
    {
    disposeStringItemList( theItem->nextString );
    disposeStringItem( theItem );
    }
  }


STRINGLIST *newStringList( void )
  {
  STRINGLIST *tempStringList;

  tempStringList = malloc( sizeof( STRINGLIST ) );
  if( tempStringList == NULL )
    {
    fprintf( stderr, "Memory allocation error in newStringList()\n" );
    return NULL;
    }

  tempStringList->firstString = NULL;
  
  return tempStringList;
  }


void disposeStringList( STRINGLIST *theStringList )
  {
  disposeStringItemList( theStringList->firstString );
  free( theStringList );
  }


ROUTINE *newRoutine( char *name, char *file )
  {
  ROUTINE *tempRoutine;

  tempRoutine = malloc( sizeof( ROUTINE ) );
  if( tempRoutine == NULL )
    {  
    fprintf( stderr, "Memory allocatin error in newRoutine()\n" );
    return NULL;
    }

  tempRoutine->name = name;
  tempRoutine->sourceFile = file;

  tempRoutine->nextRoutine = NULL;

  return tempRoutine;
  }


void disposeRoutine( ROUTINE *theRoutine )
  {
  free( theRoutine->name );
  free( theRoutine->sourceFile );
  free( theRoutine );
  }


ROUTINE *duplicateRoutine( ROUTINE *source )
  {  
  char    *tempName;
  char    *tempFile;
  ROUTINE *tempRoutine;

  tempName = malloc( strlen( source->name ) + 1 );
  if( tempName == NULL )
    {
    fprintf( stderr, "Memory allocation error in duplicateRoutine() [1]\n " );
    return NULL;
    }
  strcpy( tempName, source->name );

  tempFile = malloc( strlen( source->sourceFile ) + 1 );
  if( tempFile == NULL )
    {
    fprintf( stderr, "Memory allocation error in duplicateRoutine() [2]\n" );
    free( tempName );
    return NULL;
    }
  strcpy( tempFile, source->sourceFile );

  tempRoutine = newRoutine( tempName, tempFile );
  if( tempRoutine == NULL )
    {
    fprintf( stderr, "Memory allocation error in duplicateRoutine() [3]\n" );
    free( tempName );
    free( tempFile );
    return NULL;
    }

  return tempRoutine;
  }


INCLUDE *newInclude( char *fileName )
  {
  INCLUDE *tempInclude;
  char    *leafPtr;

  tempInclude = malloc( sizeof( INCLUDE ) );
  if( tempInclude == NULL )
    {
    fprintf( stderr, "Memory allocation error in newInclude()\n" );
    return NULL;
    }

  leafPtr = strrchr( fileName, '/' );
  if( leafPtr == NULL )
    leafPtr = fileName;
  else
    leafPtr++;

  tempInclude->name = malloc( strlen( leafPtr ) + 1 );
  if( tempInclude->name == NULL )
    {
    fprintf( stderr, "Memory allocation error in newInclude() [2]\n" );
    free( tempInclude );
    return NULL;
    }

  strcpy( tempInclude->name, leafPtr );

  tempInclude->fileName = fileName;
  tempInclude->firstCommonBlock = NULL;
  tempInclude->firstUnusedRoutine = NULL;

  tempInclude->nextInclude = NULL;
  
  return tempInclude;
  }


void disposeInclude( INCLUDE *theInclude )
  {
  free( theInclude->name );
  free( theInclude );
  }


void disposeIncludeList( INCLUDE *theInclude )
  {
  if( theInclude != NULL )
    {
    disposeIncludeList( theInclude->nextInclude );
    free( theInclude );
    }
  }


INCLUDEITEMS *newIncludeItems( void )
  {
  INCLUDEITEMS *tempIncludeItems;

  tempIncludeItems = malloc( sizeof( INCLUDEITEMS ) );
  if( tempIncludeItems == NULL )
    {
    fprintf( stderr, "Memory allocation error in newIncludeItems()\n" );
    return NULL;
    }

  tempIncludeItems->firstInclude = NULL;
  tempIncludeItems->firstInvertedRoutine = NULL;

  return tempIncludeItems;
  }


void disposeIncludeItems( INCLUDEITEMS *theIncludeItems )
  {
  disposeIncludeList( theIncludeItems->firstInclude );
  free( theIncludeItems );
  }


COMMON *newCommon( char *commonName )
  {
  COMMON *tempCommon;

  tempCommon = malloc( sizeof( COMMON ) );
  if( tempCommon == NULL )
    {
    fprintf( stderr, "Memory allocation error in newCommon()\n" );
    return NULL;
    }

  tempCommon->name = commonName;
  tempCommon->globallyUnused = False;
  tempCommon->firstUnusedRoutine = NULL;

  tempCommon->variablesSetNeverUsed = NULL;
  tempCommon->variablesUsedNeverSet = NULL;
  tempCommon->variablesNeverUsedNeverSet = NULL;

  tempCommon->nextCommonBlock = NULL;

  return tempCommon;
  }


void disposeCommon( COMMON *theCommon )
  {
  free( theCommon->name );
  free( theCommon );
  }


void disposeCommonList( COMMON *theCommon )
  {
  if( theCommon != NULL )
    {
    disposeCommonList( theCommon->nextCommonBlock );
    disposeCommon( theCommon );
    }
  }


INCLUDENAME *newIncludeName( char *name )
  {
  INCLUDENAME *tempIncludeName;
  
  tempIncludeName = malloc( sizeof( INCLUDENAME ) );
  if( tempIncludeName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newIncludeName() [1]\n" );
    return NULL;
    }

  tempIncludeName->name = malloc( strlen( name ) + 1 );
  if( tempIncludeName->name == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newIncludeName() [2]\n" );
    free( tempIncludeName );
    return NULL;
    }
  
  strcpy( tempIncludeName->name, name );
  
  tempIncludeName->nextInclude = NULL;

  return tempIncludeName;
  }


void disposeIncludeName( INCLUDENAME *theIncludeName )
  {
  free( theIncludeName->name );
  free( theIncludeName );
  }


void disposeIncludeNameList( INCLUDENAME *theIncludeName )
  {
  if( theIncludeName != NULL )
    {
    disposeIncludeNameList( theIncludeName->nextInclude );
    disposeIncludeName( theIncludeName );
    }
  }


ROUTINENAME *newRoutineName( char *name, char *file )
  {
  ROUTINENAME *tempRoutineName;

  tempRoutineName = malloc( sizeof( ROUTINENAME ) );
  if( tempRoutineName == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newRoutineName() [1]\n" );
    return NULL;
    }

  tempRoutineName->name = malloc( strlen( name ) + 1 );
  if( tempRoutineName->name == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newRoutineName() [2]\n" );
    free( tempRoutineName );
    return NULL;
    }

  strcpy( tempRoutineName->name, name );

  tempRoutineName->file = malloc( strlen( file ) + 1 );
  if( tempRoutineName->file == NULL )
    {
    fprintf( stderr, "Memory allocation failure in newRoutineName() [4]\n" );
    free( tempRoutineName->name );
    free( tempRoutineName );
    return NULL;
    }

  strcpy( tempRoutineName->file, file );
  
  tempRoutineName->firstInclude = NULL;
  tempRoutineName->nextRoutine = NULL;

  return tempRoutineName;
  }


void disposeRoutineName( ROUTINENAME *theRoutineName )
  {
  disposeIncludeNameList( theRoutineName->firstInclude );
  free( theRoutineName->name ); 
  free( theRoutineName );
  }


void disposeRoutineNameList( ROUTINENAME *theRoutineName )
  {
  if( theRoutineName != NULL )
    {
    disposeRoutineNameList( theRoutineName->nextRoutine );
    disposeRoutineName( theRoutineName );
    }
  }



/* -- Object Tests ------------------------------------------------------------
*/

BOOLEAN routineEquivilence( ROUTINE *firstRoutine, ROUTINE *secondRoutine )
  {
  if( strcmp( firstRoutine->name, secondRoutine->name ) == 0 &&
      strcmp( firstRoutine->sourceFile, secondRoutine->sourceFile ) == 0 )
    return True;
  else
    return False;
  }



/* -- List Functions ----------------------------------------------------------
*/

COMMON *lookupCommonBlock( INCLUDEITEMS *includeList, char *targetName )
  {
  INCLUDE *include;
  COMMON  *commonBlock;

  /*
   *  This badly needs rewriting as a hash table lookup
   */
  include = includeList->firstInclude;
  while( include != NULL )
    {
    commonBlock = include->firstCommonBlock;

    while( commonBlock != NULL )
      {
      if( strcmp( commonBlock->name, targetName ) == 0 )
        return commonBlock;

      commonBlock = commonBlock->nextCommonBlock;
      }

    include = include->nextInclude;
    }

  return NULL;
  } 



/* -- Queries -----------------------------------------------------------------
*/

void addRoutineToIncludeUnusedList( ROUTINE *routine, INCLUDE *includeFile )
  {
  ROUTINE *item;

  if( includeFile->firstUnusedRoutine == NULL )
    {
    includeFile->firstUnusedRoutine = duplicateRoutine( routine );
    includeFile->firstUnusedRoutine->nextRoutine = NULL;
    return;
    }
 
  /* Messy while loop and initial test condition */
  item = includeFile->firstUnusedRoutine;

  while( item->nextRoutine != NULL )
    {
    if( routineEquivilence( routine, item ) == True )
      return;
  
    item = item->nextRoutine;
    }

  if( routineEquivilence( routine, item ) == True )
    return;

  item->nextRoutine = duplicateRoutine( routine );
  item->nextRoutine->nextRoutine = NULL;
  }


BOOLEAN routineDoesntUseCommonBlock( ROUTINE *routine, COMMON *commonBlock )
  {
  ROUTINE *testRoutine;

  if( commonBlock->globallyUnused == True )
    return True;

  testRoutine = commonBlock->firstUnusedRoutine;

  while( testRoutine != NULL )
    {
    if( routineEquivilence( routine, testRoutine ) == True )
      return True;

    testRoutine = testRoutine->nextRoutine;
    }

  /* Failed to find the common block */
  return False;
  }


BOOLEAN routineDoesntUseInclude( ROUTINE *routine, INCLUDE *includeFile )
  {
  COMMON  *commonBlock;
  BOOLEAN  doesntUse;

  doesntUse = False;

  commonBlock = includeFile->firstCommonBlock;
  
  while( commonBlock != NULL )
    {
    if( routineDoesntUseCommonBlock( routine, commonBlock ) == False )
      {
      return False;
      }

    commonBlock = commonBlock->nextCommonBlock;
    }

  return True;
  }


void generateUnusedRoutinesInInclude( INCLUDE *includeFile )
  {
  COMMON  *commonBlock;
  ROUTINE *routine;

  commonBlock = includeFile->firstCommonBlock;
  
  while( commonBlock != NULL )
    {
    routine = commonBlock->firstUnusedRoutine;

    while( routine != NULL )
      {  
      if( routineDoesntUseInclude( routine, includeFile ) == True )
        addRoutineToIncludeUnusedList( routine, includeFile );

      routine = routine->nextRoutine;
      }

    commonBlock = commonBlock->nextCommonBlock;
    }
  }



/* -- Common Block Parsing Code -----------------------------------------------
*/

void parseIncludeFile( INCLUDE *include )
  {
  char       *fullIncludeName;
  FILESTATE  *file;
  PARSESTATE *state;
  char       *line;
  char       *commonName;
  COMMON     *commonBlock;

#if 0
  fullIncludeName = malloc( strlen( include->name ) + strlen( commonFilePath )
    + 1 /* terminator */ );
  strcpy( fullIncludeName, commonFilePath );
  strcat( fullIncludeName, include->name );
#endif

  file = newFileState( include->fileName );
  if( file == NULL )
    {
    fprintf( stderr, "newFileState() returned NULL in parseIncludeFile()\n" );
    printf( "The common file '%s' couldn't be opened\n", include->name );
    return;
    }

  line = getLine( file );

  while( line != NULL )
    {
    state = newParseState( line );

    if( prefix( state, "      " ) == True )
      {
      skipWord( state, "      " );
      
      if( prefix( state, "COMMON" ) == True )
        {
        skipWord( state, "COMMON" );
 
        skipUntil( state, "/" );
        skip( state, "/" );

        commonName = extract( state, "/" );

        commonBlock = newCommon( commonName );
        commonBlock->nextCommonBlock = include->firstCommonBlock;
        include->firstCommonBlock = commonBlock;
        }
      }
    
    disposeParseState( state );

    discardLine( file );
    line = getLine( file );
    }

  disposeFileState( file );
  }


void parseIncludeList( INCLUDEITEMS *includeList, char *commonListName )
  {
  FILESTATE  *file;
  PARSESTATE *state;
  char       *line;
  char       *includeName;
  INCLUDE    *include;

  file = newFileState( commonListName );
  if( file == NULL )
    {
    fprintf( stderr, "newFileState() return NULL in parseIncludeFiles()\n" );
    return;
    }

  line = getLine( file );
  while( line != NULL )
    {
    state = newParseState( line );

    skip( state, WhiteSpace );
    if( empty( state ) == False )
      {
      includeName = extract( state, WhiteSpace );

      include = newInclude( includeName );
      if( include == NULL )
        {
        fprintf( stderr, "newInclude() returned NULL in parseIncludeList()\n" );
        return;
        }

      parseIncludeFile( include );

      include->nextInclude = includeList->firstInclude;
      includeList->firstInclude = include;
      }
 
    disposeParseState( state );
    discardLine( file );
    line = getLine( file );
    }

  disposeFileState( file );
  }



/* -- Fortran Check Output Parsing --------------------------------------------
*/

void handleModuleList( COMMON *commonBlock, FILESTATE *file )
  {
  PARSESTATE *state;
  char       *routineName;
  char       *moduleName;
  char       *buffer;
  ROUTINE    *unusedRoutine;

  buffer = getLine( file );

  state = newParseState( buffer );
  if( state == NULL )
    {
    fprintf( stderr, "newParseState() return NULL in handleModuleList()\n" );
    discardLine( file );
    return;
    }

  skipWord( state, "   in module" );
  skip( state, WhiteSpace );
  routineName = extract( state, WhiteSpace );
  skipUntil( state, "]" );
  skip( state, "]" );
  moduleName = extract( state, ";""" );

  unusedRoutine = newRoutine( routineName, moduleName );
  unusedRoutine->nextRoutine = commonBlock->firstUnusedRoutine;
  commonBlock->firstUnusedRoutine = unusedRoutine;

  disposeParseState( state );

  discardLine( file );
  }


/*
 *  This routine takes too much understanding.
 */
void extractVariableStringList( STRINGLIST *stringList, FILESTATE *file, 
  char *prefixString )
  {
  PARSESTATE *state;
  char       *buffer;
  BOOLEAN     endOfList;
  char       *variable;
  STRINGITEM *stringItem;

  buffer = getLine( file );

  state = newParseState( buffer );
  if( state == NULL )
    {
    fprintf( stderr, "newParseState() returned NULL in "
      "extractVariableStringList()\n" );
    discardLine( file );
    return;
    }

  skipWord( state, prefixString );
  skip( state, WhiteSpace );

  endOfList = False;

  while( endOfList == False )
    {
    while( empty( state ) != True )
      {
      variable = extract( state, WhiteSpace );
      
      stringItem = newStringItem( variable );
      stringItem->nextString = stringList->firstString;
      stringList->firstString = stringItem;

      skip( state, WhiteSpace );
      }

    discardLine( file );
    buffer = getLine( file );
    state = newParseState( buffer );
    if( state == NULL )
      {
      fprintf( stderr, "newParseState() returned NULL in "
        "extractVariableStringList()\n" );
      discardLine( file );
      return;
      }

    if( prefix( state, " " ) != True || prefix( state, "  " ) == True )
      {
      endOfList = True;
      }
    }

  disposeParseState( state );
  }


void handleSetNeverUsed( COMMON *commonBlock, FILESTATE *file )
  {
  STRINGLIST *stringList;

  stringList = newStringList();
  if( stringList == NULL )
    {
    fprintf( stderr, "newStringList() returned NULL in "
      "handleSetNeverUsed()\n" );
    return;
    }

  extractVariableStringList( stringList, file, SetNeverUsed );

  commonBlock->variablesSetNeverUsed = stringList;
  }


void handleUsedNeverSet( COMMON *commonBlock, FILESTATE *file )
  {
  STRINGLIST *stringList;

  stringList = newStringList();
  if( stringList == NULL )
    {
    fprintf( stderr, "newStringList() returned NULL in "
      "handleSetNeverUsed()\n" );
    return;
    }

  extractVariableStringList( stringList, file, UsedNeverSet );

  commonBlock->variablesUsedNeverSet = stringList;
  }


void handleNeverUsedNeverSet( COMMON *commonBlock, FILESTATE *file )
  {
  STRINGLIST *stringList;
  
  stringList = newStringList();
  if( stringList == NULL )
    {
    fprintf( stderr, "newStringList() returned NULL in "
      "handleSetNeverUsed()\n" );
    return;
    }
  
  extractVariableStringList( stringList, file, NeverUsedNeverSet );
  
  commonBlock->variablesNeverUsedNeverSet = stringList;
  }


void parseFortranCheckCommonBlock( COMMON *commonBlock, FILESTATE *file )
  {
  PARSESTATE *state;
  char       *buffer;
  BOOLEAN     withinBlockMessages;
  BOOLEAN     moduleListFound;
  BOOLEAN     lineNotHandled;

  moduleListFound = False;
  withinBlockMessages = True;

  buffer = getLine( file );
  state = newParseState( buffer );
  skipUntil( state, ":" );
  skip( state, ": \t" );
  if( prefix( state, "unused" ) != True )
    {
    moduleListFound = True;  /* Hack - There is no list to find */
    }
  discardLine( file );

  while( withinBlockMessages == True )
    {
    buffer = getLine( file );
    if( buffer == NULL ) 
      {
      fprintf( stderr, "Unexpected end-of-file in "
        "parseFortranCheckCommonBlock()\n" );
      return;
      }

    state = newParseState( buffer );

    if( prefix( state, " " ) != True )
      {
      withinBlockMessages = False;
      if( moduleListFound == False )
        {
        commonBlock->globallyUnused = True;
        }
      }
    else
      {
      /* Test for progress */
      lineNotHandled = True; 

      if( prefix( state, SetNeverUsed ) == True )
        {
        lineNotHandled = False;
        handleSetNeverUsed( commonBlock, file );
        }

      if( prefix( state, UsedNeverSet ) == True )
        {
        lineNotHandled = False;
        handleUsedNeverSet( commonBlock, file );
        }

      if( prefix( state, NeverUsedNeverSet ) == True )
        {
        lineNotHandled = False;
        handleNeverUsedNeverSet( commonBlock, file );
        }

      if( prefix( state, "   in module" ) == True )
        {
        moduleListFound = True;
        lineNotHandled = False;
        handleModuleList( commonBlock, file );
        }

      /* Report if no progress */
      if( lineNotHandled == True )
        {
        fprintf( stderr, "line not handled in parseFortranCheckCommonBlock():"
          "\n> %s", buffer );
        discardLine( file );
        }
      }

    disposeParseState( state );
    }
  }


void parseFortranCheckOutput( INCLUDEITEMS *includeList, char *fileName )
  {
  FILESTATE  *file;
  PARSESTATE *state;
  char       *line;
  char       *commonName;
  COMMON     *commonBlock;

  file = newFileState( fileName );
  if( file == NULL )
    {
    fprintf( stderr, "newFileState() return NULL in "
      " parseFortranCheckOutput()\n");
    return;
    }

  line = getLine( file );
  while( line != NULL )
    {
    state = newParseState( line );

    if( prefix( state, "Common block " ) == True )
      {
      skipWord( state, "Common block " );
      commonName = extract( state, ":" );

      if( strcmp(commonName,"LEAD00I") && strcmp(commonName,"LEAD00R") && strcmp(commonName,"SOUR00I") && strcmp(commonName,"SOUR00R") )
      {
        commonBlock = lookupCommonBlock( includeList, commonName );
      }
      else
      {
        commonBlock = NULL;
      }
      if( commonBlock != NULL )
        {
        parseFortranCheckCommonBlock( commonBlock, file );
        }
      else
        {
        printf( "Common block %s has errors, but is not in the "
          "include lists - Please check\n",commonName );
        discardLine( file );
        }

      free( commonName ); /* Damn */
      }
    else
      {
      discardLine( file );
      }

    disposeParseState( state );

    line = getLine( file );
    }

  disposeFileState( file );
  }



/* -- HTML Output ------------------------------------------------------------
*/ 

#define Header \
  "<HTML>" \
  "<HEAD>" \
  "<TITLE>Common Block Errors</TITLE>" \
  "</HEAD>" \
  "<BODY>" \
  "<H1>Common Block Errors</H1>" \
  "<HR>" \
  "The following errors were detected by Fortran Check in the CMISS " \
  "source. Please check the routines and variables that you are reponsable " \
  "for and correct as much as possible</P>" \
  "<P>The output is broken into two seperate areas - <A HREF=\"#includes\">" \
  "Include file warnings</A> and <A HREF=\"#routines\">Routine based "\
  "warnings</A>." \
  "%s" 
  
#define Footer \
  "<HR>"

#define IncludeHeader \
  "<H3>Include <A HREF=\"" \
  CMISS_WWW_URL \
  "scripts/common-viewer.cgi" \
  "?common-file=%s\">%s</A></H3>"

#define IncludeFooter \
  "\n"


BOOLEAN variablesSetNeverUsedError( INCLUDE *include )
  {
  COMMON *commonBlock;
    
  commonBlock = include->firstCommonBlock;
  while( commonBlock != NULL )
    {
    if( commonBlock->variablesSetNeverUsed != NULL )
      return True;

    commonBlock = commonBlock->nextCommonBlock;
    }

  return False;
  }


BOOLEAN variablesUsedNeverSetError( INCLUDE *include )
  {
  COMMON *commonBlock;

  commonBlock = include->firstCommonBlock;
  while( commonBlock != NULL )
    {
    if( commonBlock->variablesUsedNeverSet != NULL )
      return True;

    commonBlock = commonBlock->nextCommonBlock;
    }

  return False;
  }


BOOLEAN variablesNeverUsedNeverSetError( INCLUDE *include )
  {
  COMMON *commonBlock;

  commonBlock = include->firstCommonBlock;
  while( commonBlock != NULL )
    {
    if( commonBlock->variablesNeverUsedNeverSet != NULL )
      return True;

    commonBlock = commonBlock->nextCommonBlock;
    }

  return False;
  }


BOOLEAN globallyUnusedError( INCLUDE *include )
  {
  COMMON *commonBlock;

  commonBlock = include->firstCommonBlock;
  while( commonBlock != NULL )
    {
    if( commonBlock->globallyUnused == True )
      return True;

    commonBlock = commonBlock->nextCommonBlock;
    }

  return False;
  }


BOOLEAN unnessessaryRoutineError( INCLUDE *include )
  {
  if( include->firstUnusedRoutine != NULL )
    return True;
  else
    return False;
  }



/*
 *  Horrible is the only way to describe this and the preceding lines.
 */
BOOLEAN includeErrors( INCLUDE *include )
  {
  if( variablesSetNeverUsedError( include ) == True )
    return True;

  if( variablesUsedNeverSetError( include ) == True )
    return True;

  if( variablesNeverUsedNeverSetError( include ) == True )
    return True;

  if( globallyUnusedError( include ) == True )
    return True;

  return False;
  }


void outputStringList( STRINGLIST *stringList, FILE *output )
  {
  STRINGITEM *stringItem;

  stringItem = stringList->firstString;
  while( stringItem != NULL )
    {
    fprintf( output, "%s ", stringItem->string );

    stringItem = stringItem->nextString;
    }
  }


void outputIncludeErrors( INCLUDE *include, FILE *output )
  {
  COMMON  *commonBlock;
  ROUTINE *routine;

  if( includeErrors( include ) == True )
    {
    fprintf( output, IncludeHeader, include->name, include->name );

    fprintf( output, "<UL>" );

    if( globallyUnusedError( include ) == True )
      {
      fprintf( output, "<LI>The following Common Blocks are globally unused in CMISS: <B>" );
      commonBlock = include->firstCommonBlock;
      while( commonBlock != NULL )
        {
        if( commonBlock->globallyUnused == True )
          fprintf( output, "%s ", commonBlock->name );

        commonBlock = commonBlock->nextCommonBlock;
        } 
      fprintf( output, "</B>" );
      }

    if( variablesSetNeverUsedError( include ) == True )
      {
      commonBlock = include->firstCommonBlock;
      while( commonBlock != NULL )
        {
        if( commonBlock->variablesSetNeverUsed != NULL )
          {
          fprintf( output, "<LI>" );
          fprintf( output, "The following variables are set but never used in Common Block %s: ", commonBlock->name );
          fprintf( output, "<B>" );
          outputStringList( commonBlock->variablesSetNeverUsed, output );
          fprintf( output, "</B>" );
          }
        
        commonBlock = commonBlock->nextCommonBlock;
        }
      }

    if( variablesUsedNeverSetError( include ) == True )
      {
      commonBlock = include->firstCommonBlock;
      while( commonBlock != NULL )
        {
        if( commonBlock->variablesUsedNeverSet != NULL )
          {
          fprintf( output, "<LI>" );
          fprintf( output, "The following variables are used but never set in Common Block %s: ", commonBlock->name );
          fprintf( output, "<B>" );
          outputStringList( commonBlock->variablesUsedNeverSet, output );
          fprintf( output, "</B>" );
          }
        
        commonBlock = commonBlock->nextCommonBlock;
        }
      }

    if( variablesNeverUsedNeverSetError( include ) == True )
      {
      commonBlock = include->firstCommonBlock;
      while( commonBlock != NULL )
        {
        if( commonBlock->variablesNeverUsedNeverSet != NULL )
          {
          fprintf( output, "<LI>" );
          fprintf( output, "The following variables are never set and never used in Common Block %s: ", commonBlock->name );
          fprintf( output, "<B>" );
          outputStringList( commonBlock->variablesNeverUsedNeverSet, output );
          fprintf( output, "</B>" );
          }
        
        commonBlock = commonBlock->nextCommonBlock;
        }
      }

    fprintf( output, "</UL>" );

    fprintf( output, IncludeFooter );
    }
  }


#define RoutineHeader \
  "<H3>Routine <A HREF=\"" \
  CMISS_WWW_URL \
  "scripts/routine-browser.cgi" \
  "?routine=%s\">%s</A> (%s)</H3>" \
  "<UL><LI>The following includes files are unnessary: "

#define RoutineFooter \
  "</UL>\n"


void printRoutineNames( FILE *output, ROUTINENAME *routine )
  {
  BOOLEAN      firstTime;
  INCLUDENAME *include;

  firstTime = True;
  include = routine->firstInclude;
  while( include != NULL )
    {
    if( firstTime == False )
      fprintf( output, ", " );
    else
      firstTime = False;
   
    fprintf( output, "<A HREF=\""
      CMISS_WWW_URL
      "scripts/common-viewer.cgi"
      "?common-file=%s\">%s</A>", include->name, include->name );
 
    include = include->nextInclude;
    }
  }


void outputRoutineErrors( ROUTINENAME *routine, FILE *output )
  {
  fprintf( output, RoutineHeader, routine->name, routine->name, 
    routine->file );

  printRoutineNames( output, routine );

  fprintf( output, RoutineFooter );
  }


void outputErrors( INCLUDEITEMS *includeList, FILE *output )
  {
  INCLUDE     *include;
  ROUTINENAME *routine;
  char         lastUpdate[256];
  time_t       systemTime;
  struct tm   *localTime;

  systemTime = time( NULL );
  localTime = localtime( &systemTime );
  strftime( lastUpdate, 255, "<P><I>This information was last updated at "
    "%H:%M on %A (%d/%m/%y)</I></P>\n", localTime );

  fprintf( output, Header, lastUpdate );

  fprintf( output, "<HR><H2><A NAME=\"includes\">Include Warnings</A></H2>"
    "<I>The following warnings list various problems with both include "
    "files and with their common blocks.</I><P>" );

  include = includeList->firstInclude;
  while( include != NULL )
    {
    outputIncludeErrors( include, output );

    include = include->nextInclude;
    }

  fprintf( output, "<HR><H2><A NAME=\"routines\">Routine Warnings</A></H2>"
    "<I>The following routines include the listed common files "
    "unnessessarily. Please remove any such include file from your "
    "routines</I><P>" );

  routine = includeList->firstInvertedRoutine;
  while( routine != NULL )
    {
    outputRoutineErrors( routine, output );

    routine = routine->nextRoutine;
    }

  fprintf( output, "<HR>" );
  }


void printUnusedList( INCLUDE *include )
  {
  ROUTINE *routine;

  routine = include->firstUnusedRoutine;
  while( routine != NULL )
    {
    printf( "routine %s does not require include file %s\n", 
      routine->name, include->name );

    routine = routine->nextRoutine;
    }
  }


void addInvertedInclude( ROUTINENAME *routine, char *includeName )
  {
  INCLUDENAME *include;
  INCLUDENAME *tempInclude;

  include = routine->firstInclude;
  while( include != NULL )
    {
    if( strcmp( include->name, includeName ) == 0 )
      return;

    include = include->nextInclude;
    }

  /* The named include wasn't found so create one */
  tempInclude = newIncludeName( includeName );
  tempInclude->nextInclude = routine->firstInclude;
  routine->firstInclude = tempInclude;
  }


void addInvertedRoutineFilePair( INCLUDEITEMS *includeList, char *includeName,
  char *routineName, char *routineFile )
  {
  ROUTINENAME *routine;
  ROUTINENAME *tempRoutine;

  routine = includeList->firstInvertedRoutine;
  while( routine != NULL )
    {
    if( strcmp( routine->name, routineName ) == 0 )
      {
      addInvertedInclude( routine, includeName );
      return;
      }

    routine = routine->nextRoutine;
    }

  /* The named routine wasn't found - create a new one */
  tempRoutine = newRoutineName( routineName, routineFile );
  tempRoutine->nextRoutine = includeList->firstInvertedRoutine;
  includeList->firstInvertedRoutine = tempRoutine;
  addInvertedInclude( includeList->firstInvertedRoutine, includeName );
  }


void invertUnusedRoutineList( INCLUDEITEMS *includeList )
  {
  INCLUDE *include;
  ROUTINE *routine;

  include = includeList->firstInclude;
  while( include != NULL )
    {
    routine = include->firstUnusedRoutine;
    while( routine != NULL )
      {
      addInvertedRoutineFilePair( includeList, include->name, routine->name, 
        routine->sourceFile );

      routine = routine->nextRoutine;
      }
    include = include->nextInclude;
    }
  }


void printInvertedUnusedRoutineList( INCLUDEITEMS *includeList )
  {
  INCLUDENAME *include;
  ROUTINENAME *routine;
    
  routine = includeList->firstInvertedRoutine;
  while( routine != NULL )
    {
    printf( "%s\n", routine->name ); 
 
    include = routine->firstInclude;
    while( include != NULL )
      {
      printf( "  %s\n", include->name );

      include = include->nextInclude;
      }
    routine = routine->nextRoutine;
    }
  }


void generateUnusedRoutines( INCLUDEITEMS *includeList )
  {
  INCLUDE *include;

  include = includeList->firstInclude;
  while( include != NULL )
    {
    generateUnusedRoutinesInInclude( include );

    include = include->nextInclude;
    }

  invertUnusedRoutineList( includeList );
  }


int lexCompareIncludeFiles( INCLUDE **firstInclude, INCLUDE **secondInclude )
  {
  return strcmp( (*firstInclude)->name, (*secondInclude)->name );
  }


void sortIncludeFiles( INCLUDEITEMS *includeList )
  {
  INCLUDE **includeArray;
  INCLUDE  *include;
  int       numElements;
  int       i;

  numElements = 0;
  include = includeList->firstInclude;
  while( include != NULL )
    {
    numElements++;
    include = include->nextInclude;
    }
  
  if( numElements < 2 )
    return;

  includeArray = malloc( numElements * sizeof( INCLUDE * ) );
  if( includeArray == NULL )
    {
    fprintf( stderr, "Memory allocation failure in sortIncludeFiles()\n" );
    return;
    }

  i = 0;
  include = includeList->firstInclude;
  while( include != NULL )
    {
    includeArray[ i++ ] = include;
    include = include->nextInclude;
    }

  qsort( includeArray, numElements, sizeof( INCLUDE * ), 
    (int (*)(const void *, const void *)) lexCompareIncludeFiles );


  includeList->firstInclude = includeArray[0];
  for( i = 0; i < numElements - 1; i++ )
    {
    includeArray[i]->nextInclude = includeArray[i+1];
    }
  includeArray[numElements - 1]->nextInclude = NULL;

  free( includeArray );  
  }


int lexCompareRoutineNames( ROUTINENAME **firstRoutine,
  ROUTINENAME **secondRoutine )
  {
  return strcmp( (*firstRoutine)->name, (*secondRoutine)->name );
  }

 
void sortInvertedFunctionList( INCLUDEITEMS *includeList )
  {
  ROUTINENAME  *routine;
  ROUTINENAME **routineArray;
  int           numRoutines;
  int           i;

  numRoutines = 0;
  routine = includeList->firstInvertedRoutine;
  while( routine != NULL )
    {
    numRoutines++;
    routine = routine->nextRoutine;
    }

  if( numRoutines < 2 )
    return;

  routineArray = malloc( numRoutines * sizeof( ROUTINENAME * ) );
  if( routineArray == NULL )
    {
    fprintf( stderr, "Memory allocation failure in "
      "sortInvertedFunctionList()\n" );
    return;
    }

  i = 0;
  routine = includeList->firstInvertedRoutine;
  while( routine != NULL )
    {
    routineArray[ i++ ] = routine;
    routine = routine->nextRoutine;
    }

  qsort( routineArray, numRoutines, sizeof( ROUTINENAME * ), 
    (int (*)(const void *, const void *)) lexCompareRoutineNames );

  includeList->firstInvertedRoutine = routineArray[0];
  for( i = 0; i < numRoutines - 1; i++ )
    {
    routineArray[i]->nextRoutine = routineArray[i+1];
    }
  routineArray[ numRoutines - 1 ]->nextRoutine = NULL;

  free( routineArray );
  }


void unnessessaryAttentionToDetailIsTheBugbearOfLittleMinds(
  INCLUDEITEMS *includeList )
  {
  sortInvertedFunctionList( includeList );
  sortIncludeFiles( includeList );
  }


void parseFortranCheck( char *fortranCheckFileName, char *commonListName,
  char *outputFileName )
  {
  INCLUDEITEMS *includeList;
  FILE         *outputFile;

  includeList = newIncludeItems();

  parseIncludeList( includeList, commonListName );
  parseFortranCheckOutput( includeList, fortranCheckFileName );
  generateUnusedRoutines( includeList );

  unnessessaryAttentionToDetailIsTheBugbearOfLittleMinds( includeList );

  outputFile = fopen( outputFileName, "wt" );
  if( outputFile == NULL )
    {
    fprintf( stderr, "fopen() failed in parseFortranCheck()\n" );
    disposeIncludeItems( includeList );
    return;
    }

  outputErrors( includeList, outputFile );

  disposeIncludeItems( includeList );
  }
