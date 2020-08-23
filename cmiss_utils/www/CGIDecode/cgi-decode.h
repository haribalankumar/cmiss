/*
 *  CGI-Decode.h
 *
 *    Public interface to the CGI decoding routines
 */
#ifndef _cgi_decode_h
#define _cgi_decode_h


/* -- Public Constants --------------------------------------------------------
*/

#define Unknown 0    /* queryType value for invalid method */
#define Get     1    /* queryType value for GET method     */
#define Post    2    /* queryType value for POST method    */



/* -- Public CGI Data Structures ----------------------------------------------
*/

#if !defined( BOOLEAN )
#define BOOLEAN int
#define True 1
#define False 0
#endif

typedef struct CGIENTRY CGIENTRY;  /* Opaque data, needed for below */

typedef struct
  {
  int       queryType;        /* POST or GET */
  float     cgiVersion;       /* e.g. CGI/1.1 == 1.1 */
  float     httpVersion;      /* e.g. HTTP/1.0 == 1.0 */
  char     *contentType;      /* MIME Content-type: == char * */
  int       contentLength;    /* MIME Content-length: == int */
  int       argc;             /* Command line argument count */
  char     *argv;             /* Command line arguments */
  
  /* Internal details - Oh for data hiding... */
  CGIENTRY *firstEntry;
  char     *queryString;
  char     *returnString;
  }
CGI;



/* -- Public Constructors/Destructors -----------------------------------------
*/

CGI *newCGI( void );
void disposeCGI( CGI *theCGI );



/* -- Public CGI Method Calls -------------------------------------------------
*/

BOOLEAN nameExists( CGI *cgi, char *name );
char *lookupString( CGI *cgi, char *name );
int lookupNumber( CGI *cgi, char *name );
CGI *getCGIEnvironment( int argc, char *argv[] );

#endif
