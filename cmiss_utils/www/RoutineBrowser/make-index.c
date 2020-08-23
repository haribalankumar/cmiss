/*
 *  Make Index.c
 *
 *    Makes the "routine-browser-list" file for the routine-browser form
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: FILE, fopen(), fgets() etc */

#include <stdlib.h>
    /* For: NULL, malloc(), free() etc */

#include <string.h>
    /* For: strlen(), etc */

#include <ctype.h>
    /* For: toupper() */

#include "wwwpaths.h"
    /* For CMISS_WWW_ROOT */

#include "simple-parser.h"
    /* For: PARSESTATE, newParseState(), prefix(), skip(), extract() etc */



/* -- Module Datatypes --------------------------------------------------------
*/

typedef struct
  {
  char *data;
  size_t startByte;
  size_t endByte;
  long  lineNumber;
  }
LINE;



/* -- Program Constants -------------------------------------------------------
*/

char *DefaultModuleFileList = "../MasterLists/module-list";
char *DefaultOutputFileName = CMISS_WWW_ROOT CM_URLPATH "/routines/routine-index";



/* -- Module Variables --------------------------------------------------------
*/

static LINE *line = NULL;
FILE        *inputFile;



/* -- Private Module Methods --------------------------------------------------
*/

char *getNextModule( FILE *moduleList )
  {
  static char theLine[256];

  fgets( theLine, 256, moduleList );

  if( feof( moduleList ) )
    return( NULL );

  if( theLine[strlen( theLine ) - 1] == '\n' )
    theLine[strlen( theLine ) - 1] = '\0';

  return theLine;
  }


void openModuleFile( char *moduleName )
  {
  inputFile = fopen( moduleName, "rt" );
  if( inputFile == NULL )
    fprintf( stderr, "File %s failed to open in openModuleFile()\n",
      moduleName );
  }


void closeModuleFile( void )
  {
  line = NULL;
  
  if( inputFile != NULL )
    fclose( inputFile );
  }


LINE *getNextLine( void )
  {
  if( inputFile == NULL )
    return NULL;

  if( line == NULL )
    {
    line = malloc( sizeof( LINE ) );
    line->data = malloc( 256 );
    line->startByte = 0;
    line->endByte = 0;
    line->lineNumber = 0;
    }

  fgets( line->data, 256, inputFile );

  if( feof( inputFile ) )
    return NULL;
  
  line->startByte = line->endByte;  /* Update from previous values */
  line->endByte = line->startByte + strlen( line->data );
  line->lineNumber++;

  return line;
  }



/*
 *  This is too ugly - fix later
 */
char *moduleID( char *moduleName )
  {
  static char  id[256];
  char        *ptr;
  char        *dest;
  char        *source;
  
  ptr = strrchr( moduleName, '/' );
  if( ptr != NULL )
    ptr++;
  else
    ptr = moduleName;

  dest = id;
  source = ptr;

  while( *source != '.' )
    {
    *dest++ = *source++;
    }
  *dest = '\0';

  ptr = id;
  while( *ptr )
    {
    *ptr = toupper( *ptr );
    *ptr++;
    }
  
  return id;
  }



#define RoutinePending 0
#define EndPending 1

void processModule( char *moduleName, FILE *outputFile )
  {
  LINE       *line;
  int         state = RoutinePending;
  PARSESTATE *parseState;
  char       *routineName;
  size_t      routineStartByte;
  size_t      routineEndByte;
  long        routineStartLine;
  long        routineEndLine;

  openModuleFile( moduleName );

  line = getNextLine();
  while( line != NULL )
    {
    parseState = newParseState( line->data );

    switch( state )
      {
      case RoutinePending:
        if( prefix( parseState, "      " ) == True )
          {
          shiftChars( parseState, 6 );

          if( prefix( parseState, "SUBROUTINE" ) == True )
            {
            skipWord( parseState, "SUBROUTINE" );
            skip( parseState, " \t" );
            routineName = extract( parseState, " (" );
            routineStartByte = line->startByte;
            routineStartLine = line->lineNumber;
            state = EndPending;
            }
          else if( prefix( parseState, "BLOCK DATA" ) == True )
            {
            skipWord( parseState, "BLOCK DATA" );
            skip( parseState, " \t" );
            routineName = extract( parseState, "\n (" );
            routineStartByte = line->startByte;
            routineStartLine = line->lineNumber;
            state = EndPending;
            }
          else
            {
            if( prefix( parseState, " " ) != True )
              {
              skipUntil( parseState, " \t" );
              skip( parseState, " \t" );
              if( prefix( parseState, "FUNCTION" ) == True )
                {
                skipWord( parseState, "FUNCTION" );
                skip( parseState, " \t" );
                routineName = extract( parseState, " (" );
                routineStartByte = line->startByte;
                routineStartLine = line->lineNumber;
                state = EndPending;
                }
              }
            }
          }
        break;

      case EndPending:
        if( prefix( parseState, "      " ) == True )
          {
          shiftChars( parseState, 6 );

          if( prefix( parseState, "END" ) == True )
            {
            skipWord( parseState, "END" );
            skip( parseState, " \t\n\r" );
            if( empty( parseState ) == True )
              {
              routineEndByte = line->endByte;
              routineEndLine = line->lineNumber;
              fprintf( outputFile, "%s %s %s %lu %ld %lu %ld\n", moduleName,
                moduleID( moduleName ),
                routineName, (unsigned long)routineStartByte, routineStartLine,
                (unsigned long)routineEndByte, routineEndLine );
              state = RoutinePending;
              }
            }
          }
        break;
      }

    disposeParseState( parseState );
    line = getNextLine();
    }

  closeModuleFile();
  }


void writeHeader( FILE *outputFile )
  {
  }


void printHelp( char *programName )
  {
  printf( "Syntax: %s [-m <module-list>] [-o <output-file]\n", programName );
  }


void makeIndex( char *outputFileName, char *moduleListName )
  {
  FILE *outputFile;
  FILE *moduleFile;
  char *moduleName;

  outputFile = fopen( outputFileName, "wt" );
  if( outputFile == NULL )
    {
    fprintf( stderr, "fopen() failed in makeIndex() [1]\n" );
    exit( 1 );
    }

  writeHeader( outputFile );

  moduleFile = fopen( moduleListName, "rt" );
  if( moduleFile == NULL )
    {
    fprintf( stderr, "fopen() failed in makeIndex() [2]\n" );
    exit( 1 );
    }

  moduleName = getNextModule( moduleFile );
  while( moduleName != NULL )
    {
    processModule( moduleName, outputFile );
    moduleName = getNextModule( moduleFile );
    }

  fclose( outputFile );
  fclose( moduleFile );
  }



/* -- Program Entry Point -----------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
  int   arg;
  char *optPtr;

  char *outputFileName = DefaultOutputFileName;
  char *moduleListName = DefaultModuleFileList;
  
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
            case 'm':
              arg++;
              if( argv[ arg ] != NULL )
                {
                moduleListName = argv[ arg ];
                }
              else
                {
                fprintf( stderr, "No Parameter for -m option\n" );
                printHelp( argv[0] );
                exit( 1 );
                }
              break;

            case 'o':
              arg++;
              if( argv[ arg ] != NULL )
                {
                outputFileName = argv[ arg ];
                }
              else
                {
                fprintf( stderr, "No Parameter for -o option\n" );
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

  makeIndex( outputFileName, moduleListName );

  return 0;
  }
