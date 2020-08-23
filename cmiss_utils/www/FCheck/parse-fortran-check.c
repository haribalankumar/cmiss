/*
 *  Parse Fortran Check.c
 *
 *    The main program for the Fortran Check Parser
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stderr, etc */

#include <stdlib.h>
    /* For: exit() */

#include "parser.h"
    /* For: parseFortranCheck() */



/* -- Program Constants -------------------------------------------------------
*/

char *DefaultCommonFileList     = "../MasterLists/common-list";
char *DefaultFortranCheckOutput = CMISS_ROOT"/cmiss/source/checkcmiss.out";
char *DefaultOutputFileName     = CMISS_WWW_SERVERPATH"help/errors/common-errors.html";



/* -- Private Functions -------------------------------------------------------
*/

void printHelp( char *programName )
  {
  printf( "Syntax: %s [-c <common file list>] [-f <fortran check output] "
    "[-o <output file>]\n", programName );
  }



/* -- Program Entry Point -----------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
  int   arg;
  char *optPtr;
  
  char *commonFileList     = DefaultCommonFileList;
  char *fortranCheckOutput = DefaultFortranCheckOutput;
  char *outputFileName     = DefaultOutputFileName;
  
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
            case 'c': 
              arg++;
              if( argv[ arg ] != NULL )
                {
                commonFileList = argv[ arg ];
                }
              else
                {
                fprintf( stderr, "No Parameter for -c option\n" );
                printHelp( argv[0] );
                exit( 1 );
                }
              break;

            case 'f':
              arg++;
              if( argv[ arg ] != NULL )
                {
                fortranCheckOutput = argv[ arg ];
                }
              else
                {
                fprintf( stderr, "No Parameter for -f option\n" );
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

  parseFortranCheck( fortranCheckOutput, commonFileList, outputFileName );
  
  return 0;
  }
