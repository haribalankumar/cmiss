/*
 *  Test-CGI.c 
 *
 *    Simple test program for the CGI-Decoer library
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: printf() */

#include <stdlib.h>
    /* For: exit() */

#include "cgi-decode.h"
    /* For: CGI handling functions */



/* -- External Prototypes -----------------------------------------------------
*/

extern int putenv(char *);



/* -- Test setup --------------------------------------------------------------
*/

void setupTestA( void )
  {
  putenv( "SERVER_SOFTWARE=NCSA/1.3" );
  putenv( "SERVER_NAME=www.esc.auckland.ac.nz" );
  putenv( "GATEWAY_INTERFACE=CGI/1.1" );
  putenv( "SERVER_PROTOCOL=HTTP/1.0" );
  putenv( "SERVER_PORT=80" );
  putenv( "REQUEST_METHOD=GET" );
  putenv( "SCRIPT_NAME=/cgi-bin/edouards-test" );
  putenv( "QUERY_STRING=customer-name=Edouard&amount-transferred=8&substring-search-flag=on&case-insensitive-flag=on&secret-password=bob+nasty" );
  putenv( "REMOTE_HOST=esu6.aukuni.ac.nz" );
  putenv( "REMOTE_ADDR=130.216.5.106" );
  }


void setupTestB( void )
  {
  putenv( "SERVER_SOFTWARE=NCSA/1.3" );
  putenv( "SERVER_NAME=www.esc.auckland.ac.nz" );
  putenv( "GATEWAY_INTERFACE=CGI/1.1" );
  putenv( "SERVER_PROTOCOL=HTTP/1.0" );
  putenv( "SERVER_PORT=80" );
  putenv( "REQUEST_METHOD=GET" );
  putenv( "SCRIPT_NAME=/cgi-bin/edouards-test" );
  putenv( "QUERY_STRING=340,124" );
  putenv( "REMOTE_HOST=esu6.aukuni.ac.nz" );
  putenv( "REMOTE_ADDR=130.216.5.106" );
  }




/* -- The Program Entry Point -------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
  CGI *cgi;
  char *customerName;
  int   amountTransferred;
  char *thePassword;
  int   x, y;

  printf( "Test A (normal CGI test)\n" );
 
  setupTestA();

  cgi = getCGIEnvironment( argc, argv );
  if( cgi == NULL )
    {
    fprintf( stderr, "getCGIEnvironment() failed in main()\n" );
    exit( 1 );
    }

  customerName = lookupString( cgi, "customer-name" );
  amountTransferred = lookupNumber( cgi, "amount-transferred" );

  printf( "%s transferred %d\n", customerName, amountTransferred );

  if( nameExists( cgi, "secret-password" ) == True )
    {
    thePassword = lookupString( cgi, "secret-password" );
    printf( "The secrtet password for today is: %s\n", thePassword );
    }

  disposeCGI( cgi );

  printf( "\nTest B (map test)\n" );

  setupTestB();

  cgi = getCGIEnvironment( argc, argv );
  if( cgi == NULL )
    {
    fprintf( stderr, "getCGIEnvironment() failed in main()\n" );
    exit( 1 );
    }

  x = lookupNumber( cgi, "x" );
  y = lookupNumber( cgi, "y" );

  printf( "Map position is %d, %d\n\n", x, y );

  disposeCGI( cgi );

  return 0;
  }

