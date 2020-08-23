/*
 * List Dir.c
 *
 *    Lists the users input or output directories.
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

#include "cgi-decode.h"
    /* For: CGI handling functions */

#include "groups.h"
    /* For: GROUP, newGroup(), getTag(), etc */

#include "tags.h"
    /* For: TAG, newTag(), etc */

#include "rule-parse.h"
    /* For: parseRuleFile() */

#include "user-managment.h"
    /* For: USER, USERLIST, readUserList(), passwordCheck(), etc */

#include "get-user.h"
    /* For: getUser() */

#include "job.h"
    /* For: JOB, getJob(), printJobHeading(), etc */

#include "directory.h"
    /* For: FILEINFO, DIRECTORY, newDirectory() etc */

#include "border.h"
    /* For: getBorderHeader() */

#include "webcmiss.h"


/* -- Module Data Structures --------------------------------------------------
*/

/* Bitmask of possible file actions */
typedef int FILEACTIONS;

/* Bitmask Constants */
const int EditFileAction   = 0x0001;
const int RevertFileAction = 0x0002;
const int ViewFileAction   = 0x0004;
const int DeleteFileAction = 0x0008;
const int InfoFileAction   = 0x0010;
const int CmguiAction      = 0x0020;
const int PluginAction     = 0x0040;

#define BUFFER_SIZE 1024


/* -- Method Prototypes -------------------------------------------------------
*/

void addFileInfoToDirectory( FILEINFO *fileInfo, DIRECTORY *directory );



/* -- Test Code ---------------------------------------------------------------
*/ 

#ifdef DEBUG
static void test( char *pathName )
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
#endif

/* -- Support Code ------------------------------------------------------------
*/ 

/*
 *  dup()
 */
char *newStringCopy( const char *source )
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



/*
 *  System Specific (i.e. this only work under UNIX) 
 */
char *newLeafString( const char *path )
  {
  char       *tempString;
  const char *ptr;

  /* Step 1. find leaf part of name */
  ptr = strrchr( path, '/'  );
  if( ptr == NULL )
    ptr = path;
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



/* -- URL utilities -----------------------------------------------------------
*/ 

char *makeEditURL( USER *user, char *file )
  {
  char *tempString;
  char *name;
  char *password;

  char *EditFileProgram = EDIT_FILE;

  name = user->name;
  password = user->password;

  tempString = malloc( 
    strlen( EditFileProgram ) + 1 +       /* "Scripts/blah" + "?" */
    5 + strlen( name ) +                  /* "user=" + "fred"     */
    1 +                                   /* "&"                  */
    9 + strlen( password ) +              /* "password" + "frob"  */  
    1 +                                   /* "&"                  */
    5 + strlen( file ) +                  /* "file=" + "blem.com" */
    1 );                                  /* "\0"                 */

  if( tempString == NULL )
    {
    fprintf( stderr, "Memory allocation failure in makeEditURL()\n" );
    return NULL;
    }

  /* NB. I don't url-escape this string yet. Do so later */
  sprintf( tempString, "%s?user=%s&password=%s&file=%s", EditFileProgram, 
    name, password, file );

  return tempString;
  }


char *makeRevertURL( USER *user, char *file )
  {
  char *tempString;
  char *name;
  char *password;

  char *RevertFileProgram = REVERT_FILE;

  name = user->name;
  password = user->password;

  tempString = malloc( 
    strlen( RevertFileProgram ) + 1 +     /* "Scripts/blah" + "?" */
    5 + strlen( name ) +                  /* "user=" + "fred"     */
    1 +                                   /* "&"                  */
    9 + strlen( password ) +              /* "password" + "frob"  */  
    1 +                                   /* "&"                  */
    5 + strlen( file ) +                  /* "file=" + "blem.com" */
    1 );                                  /* "\0"                 */

  if( tempString == NULL )
    {
    fprintf( stderr, "Memory allocation failure in makeEditURL()\n" );
    return NULL;
    }

  /* NB. I don't url-escape this string yet. Do so later */
  sprintf( tempString, "%s?user=%s&password=%s&file=%s", RevertFileProgram, 
    name, password, file );

  return tempString;
  }


char *makeViewURL( USER *user, char *file, char *directory )
  {
  char *tempString;
  char *name;
  char *password;

  char *ViewFileProgram = VIEW_FILE;

  name = user->name;
  password = user->password;

  tempString = malloc( 
    strlen( ViewFileProgram ) + 1 +     /* "Scripts/blah" + "?" */
    5 + strlen( name ) +                  /* "user=" + "fred"     */
    1 +                                   /* "&"                  */
    9 + strlen( password ) +              /* "password" + "frob"  */  
    1 +                                   /* "&"                  */
   10 + strlen( directory ) +             /* "directory" + "blob" */  
    1 +                                   /* "&"                  */
    5 + strlen( file ) +                  /* "file=" + "blem.com" */
    1 +                                   /* "&"                  */
    7 + strlen( LIST_PROBLEM ) +          /* "return=" + "list-problem.cgi" */
    1 );                                  /* "\0"                 */

  if( tempString == NULL )
    {
    fprintf( stderr, "Memory allocation failure in makeViewURL()\n" );
    return NULL;
    }

  /* NB. I don't url-escape this string yet. Do so later */
  sprintf( tempString, "%s?user=%s&password=%s&directory=%s&file=%s&return=%s", 
					ViewFileProgram, name, password, directory, file, LIST_PROBLEM );

  return tempString;
  }


char *makeDeleteURL( USER *user, char *file, char *directory )
  {
  char *tempString;
  char *name;
  char *password;

  char *DeleteFileProgram = DELETE_FILE;

  name = user->name;
  password = user->password;

  tempString = malloc( 
    strlen( DeleteFileProgram ) + 1 +     /* "Scripts/blah" + "?" */
    5 + strlen( name ) +                  /* "user=" + "fred"     */
    1 +                                   /* "&"                  */
    9 + strlen( password ) +              /* "password" + "frob"  */  
    1 +                                   /* "&"                  */
   10 + strlen( directory ) +             /* "directory" + "blob" */  
    1 +                                   /* "&"                  */
    5 + strlen( file ) +                  /* "file=" + "blem.com" */
    1 );                                  /* "\0"                 */

  if( tempString == NULL )
    {
    fprintf( stderr, "Memory allocation failure in makeDeleteURL()\n" );
    return NULL;
    }

  /* NB. I don't url-escape this string yet. Do so later */
  sprintf( tempString, "%s?user=%s&password=%s&directory=%s&file=%s", 
					DeleteFileProgram, name, password, directory, file );
  return tempString;
  }

char *makeCmguiURL( USER *user, char *file )
  {
  char *tempString;
  char *name;
  char *password;

  char *DeliverFileProgram = CMGUI_PROGRAM;

  name = user->name;
  password = user->password;

  tempString = malloc( 
    strlen( DeliverFileProgram ) + 1 +    /* "Scripts/blah" + "?" */
    5 + strlen( name ) +                  /* "user=" + "fred"     */
    1 +                                   /* "&"                  */
    9 + strlen( password ) +              /* "password" + "frob"  */  
    1 +                                   /* "&"                  */
    5 + strlen( file ) +                  /* "file=" + "blem.com" */
    1 );                                  /* "\0"                 */

  if( tempString == NULL )
    {
    fprintf( stderr, "Memory allocation failure in makeCmguiURL()\n" );
    return NULL;
    }

  /* NB. I don't url-escape this string yet. Do so later */
  sprintf( tempString, "%s?user=%s&password=%s&file=%s", DeliverFileProgram, 
    name, password, file );
  return tempString;
  }

char *makeFormURL( USER *user, char *file )
  {
  char *tempString;
  char *name;

  name = user->name;

  tempString = malloc( strlen( name ) + strlen( WEBCMISS_URL ) + 35);
  if( tempString == NULL )
    {
    fprintf( stderr, "Memory allocation failure in makeFormURL()\n" );
    return NULL;
    }

  /* NB. I don't url-escape this string yet. Do so later */
  sprintf( tempString, "%speople/%s/Input/%s", WEBCMISS_URL, name, file);

  return tempString;
  }


/* -- Table Code --------------------------------------------------------------
*/

void printTableSurroundStart( FILE *output )
  {
  fprintf( output, "<TABLE WIDTH=\"100%%\" BORDER=\"0\" "
    "CELLSPACING=\"0\" CELLPADDING=\"0\">" );
  }


void printHeaderInfo( FILE *output, char *title, char *url )
  {
  fprintf( output, "<TR>" );
  fprintf( output, "<TD ALIGN=\"LEFT\" VALIGN=\"BOTTOM\">"
    "<IMG SRC=\"" IMAGE_DIR "folder.gif\" WIDTH=32 HEIGHT=32>"
    "<FONT SIZE=\"+2\">%s</FONT></TD>", title );
	if (url != NULL )
		{
    fprintf( output, "<TD ALIGN=\"RIGHT\" VALIGN=\"BOTTOM\">"
						"<FONT SIZE=\"-1\"> (<A HREF=\"%s\">View Example Tree</A>)"
						"</FONT></TD>", url );
	  }
  fprintf( output, "</TR>" );
  }


void printIndentedTableStart( FILE *output )
  {
  fprintf( output, "<TR><TD COLSPAN=\"2\">"
    "<TABLE WIDTH=\"100%%\" BORDER=\"0\" CELLSPACING=\"0\" CELLPADDING=\"0\">");
  fprintf( output, "<TR>" );
  fprintf( output, "<TD WIDTH=\"10%%\"> </TD> <TD WIDTH=\"90%%\">" );
  }


void printFiveColumnStart( FILE *output )
  {
  fprintf( output,
    "<TABLE WIDTH=\"100%%\" BORDER=\"0\" CELLSPACING=\"4\" CELLPADDING=\"0\">"
    "<TD WIDTH=\"100%%\" COLSPAN=\"5\"> <HR NOSHADE> </TD>" );
  }


void printColumnNames( FILE *output )
  {
  fprintf( output, "<TR><TH>Name<HR NOSHADE></TH><TH>Action<HR NOSHADE></TH>"
    "<TH>Modified<HR NOSHADE></TH><TH>Size<HR NOSHADE></TH>"
    "<TH>Description<HR NOSHADE></TH></TR>\n" );
  }


void printEmptyDirectory( FILE *output )
  {
  fprintf( output, "<TR><TD COLSPAN=\"5\" ALIGN=\"CENTER\">"
    "<B>The directory is empty</B></TD></TR>\n" );
  }


void printFileNameAndIcon( FILE *output, FILEINFO *file, USER *user, 
  FILEACTIONS fileActions, char *directory )
  {
  char *URL;

  if(!strstr(file->name,".form"))
  {
    if( fileActions & EditFileAction )
    URL = makeEditURL( user, file->name );
    else if( fileActions & ViewFileAction )
    URL = makeViewURL( user, file->name, directory ); 
    else
    URL = NULL;
  }
  else
  {
  URL = makeFormURL( user, file->name ); 
  }
  
  if( file->type == DirectoryType )
    fprintf( output, "<TD ALIGN=\"CENTER\">"
      "<IMG SRC=\"" IMAGE_DIR "folder.gif\" WIDTH=32 HEIGHT=32><BR><I>%s</I>"
      "</TD>", file->name );
  else
    if( URL != NULL )
      fprintf( output, "<TD ALIGN=\"CENTER\">"
        "<A HREF=\"%s\">"
        "<IMG BORDER=\"0\" SRC=\"" IMAGE_DIR "file.gif\" WIDTH=29 HEIGHT=32>"
        "<BR><I>%s</I></A></TD>",
        URL, file->name );
    else
      fprintf( output, "<TD ALIGN=\"CENTER\">"
        "<IMG SRC=\"" IMAGE_DIR "file.gif\" WIDTH=29 HEIGHT=32>"
        "<BR><I>%s</I></A></TD>", file->name );

  if( URL != NULL )
    free( URL );
  }


void printFileActions( FILE *output, FILEINFO *file, FILEACTIONS fileActions,
  USER *user, char *directory )
  {
    int i;
    char *actionURL;

  /* No file actions for a directory */
  if( file->type == DirectoryType )
    {
    fprintf( output, "<TD> </TD>" );
    return;
    }

  fprintf( output, "<TD ALIGN=\"CENTER\">" );

  if( fileActions & EditFileAction )
  { 
    actionURL = makeEditURL( user, file->name ); 
    if( actionURL != NULL ) 
      fprintf( output, "<A HREF=\"%s\">", actionURL ); 
    fprintf( output, "<IMG vspace=\"1\" BORDER=\"0\" SRC=\"" IMAGE_DIR "edit.gif\""
      " WIDTH=38 HEIGHT=12><BR>" ); 
    if( actionURL != NULL ) 
      fprintf( output, "</A>" ); 
    if( actionURL != NULL ) 
      free( actionURL ); 
  } 

  if( fileActions & RevertFileAction )
  {
    actionURL = makeRevertURL( user, file->name ); 
    if( actionURL != NULL ) 
      fprintf( output, "<A HREF=\"%s\">", actionURL ); 
    fprintf( output, "<IMG VSPACE=\"1\" BORDER=\"0\" SRC=\"" IMAGE_DIR "revert.gif\""
      " WIDTH=38 HEIGHT=12><BR>" ); 
    if( actionURL != NULL ) 
      fprintf( output, "</A>" ); 
    if( actionURL != NULL ) 
      free( actionURL ); 
  }

  if( fileActions & ViewFileAction )
  { 
    actionURL = makeViewURL( user, file->name, directory ); 
    if( actionURL != NULL ) 
      fprintf( output, "<A HREF=\"%s\">", actionURL ); 
    fprintf( output, "<IMG vspace=\"1\" BORDER=\"0\" SRC=\"" IMAGE_DIR "view.gif\""
      " WIDTH=38 HEIGHT=12><BR>" ); 
    if( actionURL != NULL ) 
      fprintf( output, "</A>" ); 
    if( actionURL != NULL )
      free( actionURL ); 
  } 

  if( fileActions & DeleteFileAction )
  { 
    actionURL = makeDeleteURL( user, file->name, directory ); 
    if( actionURL != NULL ) 
      fprintf( output, "<A HREF=\"%s\">", actionURL ); 
    fprintf( output, "<IMG VSPACE=\"1\" BORDER=\"0\" SRC=\"" IMAGE_DIR "delete.gif\""
      " WIDTH=38 HEIGHT=12><BR>" ); 
    if( actionURL != NULL ) 
      fprintf( output, "</A>" ); 
    if( actionURL != NULL ) 
      free( actionURL ); 
  } 

  if( fileActions & InfoFileAction )
    fprintf( output, "<IMG VSPACE=\"1\" BORDER=\"0\" SRC=\"" IMAGE_DIR "info.gif\""
      " WIDTH=38 HEIGHT=12><BR>" );

  if( fileActions & CmguiAction )
    {
		for(i=0;i<strlen( file->name )-3;i++)
			{
      if(!strncmp(file->name+i,".com",4))
				{
        actionURL = makeCmguiURL( user, file->name );  
        if( actionURL != NULL )
        fprintf( output, "<A HREF=\"%s\">", actionURL );
        fprintf( output, "<IMG VSPACE=\"1\" BORDER=\"0\" SRC=\"" IMAGE_DIR "cmgui.gif\""
          " WIDTH=38 HEIGHT=12><BR>" );
        if( actionURL != NULL )
        fprintf( output, "</A>" );
        if( actionURL != NULL )
        free( actionURL );
        }
      }
	  }

	fprintf( output, "</TD>" );
  }

void printFileTimeModified( FILE *output, FILEINFO *file )
  {
  time_t     timeNow;
  time_t     timeModified;
  struct tm *brokenDownTime;
  double     difference;
  char       timeBuffer[64];

  timeNow        = time( NULL );
  timeModified   = file->modificationDate;
  brokenDownTime = localtime( &timeModified );
  difference     = difftime( timeNow, timeModified );

  strftime( timeBuffer, 64, "%I:%M%p %d/%m/%y", brokenDownTime );

  fprintf( output, "<TD ALIGN=\"CENTER\"><FONT SIZE=\"-1\">" );
  fprintf( output, "%s", timeBuffer );
  fprintf( output, "<BR>" );

  fprintf( output, "(" );

  if( difference < 2.0 )
    fprintf( output, "1 second ago" );
  else if( difference < 60.0 )
    fprintf( output, "%d seconds ago", (int) difference );
  else if( difference < 120.0 )
    fprintf( output, "1 minute ago" );
  else if( difference < 3600.0 )
    fprintf( output, "%d minutes ago", (int) (difference/60.0) );
  else if( difference < 7200.0 )
    fprintf( output, "1 hour ago" );
  else if( difference < 86400.0 )
    fprintf( output, "%d hours ago", (int) (difference/3600.0) );
  else if( difference < 172800.0 )
    fprintf( output, "1 day ago" );
  else if( difference < 604800.0 )
    fprintf( output, "%d days ago", (int) (difference/86400.0) );
  else if( difference < 1209600.0 )
    fprintf( output, "1 week ago" );
  else if( difference < 31536000.0 )
    fprintf( output, "%d weeks ago", (int) (difference/604800.0) );
  else if( difference < 63072000.0 ) 
    fprintf( output, "1 year ago" );
  else
    fprintf( output, "%d years ago", (int) (difference/31536000.0) );

  fprintf( output, ")" );

  fprintf( output, "</FONT></TD>" );
  }


void printFileNeverModified( FILE *output )
  {
  fprintf( output, "<TD ALIGN=\"CENTER\"><FONT SIZE=\"-1\">" );
  fprintf( output, "Original File" );
  fprintf( output, "</FONT></TD>" );
  }


/* NB This will fail on a 16-bit machine. Fix later */
void printFileSize( FILE *output, int size )
  {
  fprintf( output, "<TD ALIGN=\"CENTER\"><FONT SIZE=\"-1\">" );

  if( size < 1024 )
    fprintf( output, "%d bytes", size );
  else if( size < 2048 )
    fprintf( output, "1 K" );
  else if( size < 1048576 )
    fprintf( output, "%d K", size/1024 );
  else if( size < 2097152 )
    fprintf( output, "1 M" );
  else if( size < 1073741824 )
    fprintf( output, "%d M", size/1048576 );
  else if( size < 2147483648 )
    fprintf( output, "1 G" );
  else 
    fprintf( output, "%d G", size/1073741824 );

  fprintf( output, "</FONT></TD>" );
  }


void printFileDescription( FILE *output, char *description )
  {
  fprintf( output, "<TD>%s</TD>", description );
  }


void printFileInfo( FILE *output, FILEINFO *file, USER *user,
  FILEACTIONS fileActions, char *directory )
  {
    if(strncmp(file->name,"exported",8)&&strncmp(file->name,"description.html",16)&&
      strncmp(file->name,"name.txt",8)&&!strstr(file->name,".cgi")
      &&!strstr(file->name,"_template")&&!strstr(file->name,"gif"))
    {
      fprintf( output, "<TR>" );
      
      printFileNameAndIcon( output, file, user, fileActions, directory );
      
      printFileActions( output, file, fileActions, user, directory );
      
      if( file->modificationDate == (time_t) -1 )
        {
        printFileNeverModified( output );
        }
      else
        {
        printFileTimeModified( output, file );
        }
      
      printFileSize( output, file->size );
      
      printFileDescription( output, file->description );
      
      fprintf( output, "</TR>\n" );
    }
  }


void printFiveColumnEnd( FILE *output )
  {
  fprintf( output, "<TR><TD COLSPAN=\"5\"><HR NOSHADE></TD></TR></TABLE>" );
  }


void printIndentedTableEnd( FILE *output )
  {
  fprintf( output, "</TD></TR></TR></TABLE>" );
  }


void printTableSurroundEnd( FILE *output )
  {
  fprintf( output, "</TABLE>" );
  }


void printChangeButton( FILE *output, USER *user )
  {
  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" CHANGE_EXAMPLE "\">" );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">",
    user->name );
  fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">", 
    user->password );
  fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"Switch Example\">" );
  fprintf( output, "</FORM>" );
  }


void printDeleteButton( FILE *output, USER *user, char *directory )
  {
  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" DELETE_DIR "\">" );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">",  user->name );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"password\" VALUE=\"%s\">", user->password );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"directory\" VALUE=\"%s\">",	directory );
  fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"Delete Files\">" );
  fprintf( output, "</FORM>" );
  }


void printPasswordButton( FILE *output, USER *user )
  {
  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" CHANGE_PASSWORD "\">" );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">",
    user->name );
  fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">", 
    user->password );
  fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"Change Password\">" );
  fprintf( output, "</FORM>" );
  }


void printHTMLDirectory( FILE *output, DIRECTORY *directory, char *directoryTitle, 
												 char *directoryName, USER *user, FILEACTIONS fileActions, 
												 BOOLEAN inputDir )
  {
  FILEINFOITEM *item;
  FILEINFO     *file;
	char *examplesURL = NULL;

	if (inputDir)
		{
		examplesURL = EXAMPLE_FILES_URL;
		}

  printTableSurroundStart( output );
	printHeaderInfo( output, directoryTitle, examplesURL );
  printIndentedTableStart( output );
  printFiveColumnStart( output );

  item = directory->fileInfoList->firstFileInfoItem;

  if( item == NULL )
    {
    printEmptyDirectory( output );
    }
  else
    {
    printColumnNames( output );

    while( item != NULL )
      {
      file = item->fileInfo;
      printFileInfo( output, file, user, fileActions, directoryName );
      item = item->nextFileInfoItem;
      }
    }

  printFiveColumnEnd( output );
  printIndentedTableEnd( output );
  printTableSurroundEnd( output );
  }


void listProblemDirectory( CGI *cgi, FILE *output, BOOLEAN inputDir )
  {
  USER      *user;
  JOB       *job;

  /* Constants */
  const FILEACTIONS InputFileActions =
    (EditFileAction | InfoFileAction | RevertFileAction);
  const FILEACTIONS OutputFileActions = 
    (ViewFileAction | DeleteFileAction | PluginAction | CmguiAction);

  user = getUser( cgi );
  if( user == NULL )
    {
    /* Not always an error, the user could simply have supplied an incorrect
       password. getUser() will report back those conditions to the user
       (one hopes...) */
    return;
    }

  job = getJob( user );
  if( job == NULL )
    {
    fprintf( stderr, "getJob() returned NULL in listProblem()\n" );
    disposeUser( user );
    return;
    }


  printJobHeading( output, user, job );
  fprintf( stdout, "<P>\n" );

	if ( inputDir )
		{
		printHTMLDirectory( output, job->inputFiles, "Input Files", "input" , user,
											 InputFileActions, inputDir );
		fprintf( output, "<P>\n" );

		fprintf( output, "<CENTER>\n" );
		fprintf( output, "<TABLE BORDER=\"0\" WIDTH=\"70%%\">\n" );
		fprintf( output, "<TR><TD ALIGN=\"LEFT\" WIDTH=\"25%%\">\n" );
		printDeleteButton( output, user, "input" );
		fprintf( output, "</TD>\n" );
		fprintf( output, "<TD ALIGN=\"LEFT\" WIDTH=\"25%%\">\n" );
		printChangeButton( output, user );
		fprintf( output, "</TR>\n</TABLE>\n" );
		fprintf( output, "</CENTER>\n" );
 
		fprintf( output, "<P>\n" );

	  }
	else 
		{
			printHTMLDirectory( output, job->outputFiles, "Output Files", "output", user, 
												 OutputFileActions, inputDir );
			fprintf( output, "<P>\n" );

			fprintf( output, "<CENTER>\n" );
			printDeleteButton( output, user, "output" );
			fprintf( output, "</CENTER>\n" );
		}
  /* Free stuff (todo) */
  }



/* -- Program Entry Point -----------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
  CGI *cgi;
	USER *user;
	BOOLEAN inputDir;
	char *dir;


  cgi = getCGIEnvironment( argc, argv );
  if( cgi == NULL )
    {
    printf( "<HR><P><B>ERROR:</B> getCGIEnvironment() failed in main()</P>"
            "<P>Email "WEBMASTER"<I>immediately</I></P>\n" );
    return 0;
    }

  user = getUser( cgi );

  if(!(dir = lookupString(cgi,"directory")))
    {
    dir = "input";
    }

	if ( strcmp(dir,"input") == 0 )
		{
		inputDir = True;
	  }
	else
		{
		inputDir = False;
		}


	borderGenerateHeader( stdout, user, False, NULL, NULL );
  
  listProblemDirectory( cgi, stdout, inputDir );

	borderGenerateFooter( stdout );

  return 0;
  }

