/*
 *  Account Details.c
 *
 *    Dialog of account options.
 */


/* -- Include Files -----------------------------------------------------------
*/

#include <stdio.h>
    /* For: fprintf(), stderr, etc */

#include <stdlib.h>
    /* For: malloc(), free(), NULL, etc */

#include <string.h>
    /* For: strlen(), strcpy(), strstr(), etc */

#include "get-user.h"
    /* For: getUser() */

#include "job.h"
    /* For: JOB, getJob(), printJobHeading(), etc */

#include "border.h"
    /* For: borderGenerateHeader etc */

#include "host.h"
    /* For: HOSTLIST, etc */

#include "resource.h"
    /* For: getDiskUsage(), fprintDiskUsage() */

#include "webcmiss.h"

/* -- Module Methods ----------------------------------------------------------
*/

void printLogButton( FILE *output, USER *user )
  {
  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" VIEW_FILE "\">" );
  fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">", user->name );
  fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">", user->password );
  fprintf( output, "<INPUT NAME=\"directory\" TYPE=\"hidden\" VALUE=\"root\">" );
  fprintf( output, "<INPUT NAME=\"file\" TYPE=\"hidden\" VALUE=\"job.log\">" );
  fprintf( output, "<INPUT NAME=\"return\" TYPE=\"hidden\" VALUE=\"" ACCOUNT_DETAILS "\">" );
  fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"View Log File\">" );
  fprintf( output, "</FORM>" );
  }

void printOKButton( FILE *output, USER *user )
  {
  fprintf( output, "<FORM METHOD=\"POST\" ACTION=\"" ACCOUNT_DETAILS "\">" );
  fprintf( output, "<INPUT NAME=\"user\" TYPE=\"hidden\" VALUE=\"%s\">", user->name );
  fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">", user->password );
  fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"OK\">" );
  fprintf( output, "</FORM>" );
  }

void printResourceUsage( FILE *output, USER *user, HOSTLIST *hostList )
  {
  JOB        *job;
	USAGELIST  *usageList;
  char       buff[1024];
	int        i, nhost;

  job = getJob( user );
  if( job == NULL )
    {
    fprintf( stderr, "getJob() returned NULL in abortJob()\n" );
    disposeUser( user );
    return;
    }

  sprintf( buff, "%s/job.log", user->homeDirectory );
	usageList = getUsageFromLog(buff, hostList);
	nhost = usageList->nhost;

  
  fprintf( output, "<H1>Resource Usage</H1><P>\n" );
  fprintf( output, "<HR NOSHADE><P>\n" );

  /*
   * Job Status
   */
  fprintf( output, "<H3>Job Status</H3><P>\n" );
  printJobHeading( output, user, job );
  fprintf( output, "<P><HR NOSHADE>\n" );

  /*
   * CPU Usage
   */
  fprintf( output, "<H3>CPU Usage</H3><P>\n" );
  fprintf( output, "<CENTER>\n" );

  fprintf( output, "<FONT SIZE=\"+1\">CPU Usage to date by Host</FONT><P>\n" );

  fprintf( output, "<TABLE BORDER>\n" );
  fprintf( output, "  <TR>\n" );
  fprintf( output, "    <TH>Hostname</TH>\n" );
  fprintf( output, "    <TH>Wall Time</TH>\n" );
  fprintf( output, "    <TH>CPU Time</TH>\n" );
  fprintf( output, "  </TR>\n" );

	for (i=0;i<nhost;i++) {
    fprintf( output, "  <TR>\n" );
    fprintf( output, "    <TH>%s</TH>\n",   usageList->usage[i].hostname );
    fprintf( output, "    <TD>%.1f</TD>\n", usageList->usage[i].wall_total );
    fprintf( output, "    <TD>%.1f</TD>\n", usageList->usage[i].cpu_total );
    fprintf( output, "  </TR>\n" );
    }
  fprintf( output, "</TABLE><P>\n" );
  
  fprintf( output, "<FONT SIZE=\"+1\">Log of last job by Host</FONT><P>\n" );
  fprintf( output, "<TABLE BORDER>\n" );

  fprintf( output, "  <TR>\n" );
  fprintf( output, "    <TH>Host</TH>\n" );
  fprintf( output, "    <TH>Start Time</TH>\n" );
  fprintf( output, "    <TH>Wall Time</TH>\n" );
  fprintf( output, "    <TH>CPU Time</TH>\n" );
  fprintf( output, "    <TH>Version</TH>\n" );
  fprintf( output, "    <TH>Threads</TH>\n" );
  fprintf( output, "    <TH>Description</TH>\n" );
  fprintf( output, "  </TR>\n" );

	for (i=0;i<nhost;i++) {
    fprintf( output, "  <TR>\n" );
    fprintf( output, "    <TH>%s</TH>\n",   usageList->usage[i].hostname );
    fprintf( output, "    <TD>%s</TD>\n",   usageList->usage[i].start_time );
    if (usageList->usage[i].start_time[0] != '\0') {
			fprintf( output, "    <TD>%.1f</TD>\n", usageList->usage[i].wall_time );
			fprintf( output, "    <TD>%.1f</TD>\n", usageList->usage[i].cpu_time );
		}
		else {
			fprintf( output, "    <TD></TD>\n" );
			fprintf( output, "    <TD></TD>\n" );
		}
    fprintf( output, "    <TD>%s</TD>\n", usageList->usage[i].version );
    fprintf( output, "    <TD>%s</TD>\n", usageList->usage[i].threads );
    fprintf( output, "    <TD>%s</TD>\n", usageList->usage[i].description );
    fprintf( output, "  </TR>\n" );
    }
  fprintf( output, "</TABLE><P>\n" );

  fprintf( output, "<P>\n" );
  printLogButton( output, user );
  fprintf( output, "<P>\n" );

  fprintf( output, "</CENTER>\n" );
  fprintf( output, "<P><HR NOSHADE><P>\n" );

  /*
   * Disk Usage
   */
  fprintf( output, "<H3>Disk Usage</H3><P>\n" );
  fprintf( output, "<CENTER>\n" );
  fprintf( output, "<FONT SIZE=\"+1\">The account is using " );
  fprintDiskUsage( output, getDiskUsage(user->homeDirectory) );
  fprintf( output, " of disk space.</FONT><P>\n" );
  
  fprintf( output, "<TABLE><TR>\n" );
  fprintf( output, "<TH>Input directory:  <TD>" );
  sprintf( buff, "%s/Input", user->homeDirectory );
  fprintDiskUsage( output, getDiskUsage(buff) );
  fprintf( output, "\n<TR>" );
  fprintf( output, "<TH>Output directory: <TD>" );
  sprintf( buff, "%s/Output", user->homeDirectory );
  fprintDiskUsage( output, getDiskUsage(buff) );
  fprintf( output, "</TABLE><P>\n" );
  fprintf( output, "</CENTER>\n" );
  fprintf( output, "<P><HR NOSHADE><P>\n" );

  fprintf( output, "<CENTER>" );
  printOKButton( output, user );
  fprintf( output, "</CENTER>" );

	freeUsageList(usageList);
  }


void printChangeForm( FILE *output, USER *user )
  {

  fprintf( output, "<H1>Change Password</H1><P>" );
  fprintf( output, "<HR NOSHADE><P>" );

  fprintf( output, "<FONT SIZE=\"+1\">Change the password for user \"%s\"</FONT><P>\n", user->name );

	fprintf( output, "<FORM ACTION=\"" CHANGE_PASSWORD "\" METHOD=\"POST\">\n" );
	fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">\n", user->name );
	fprintf( output, "<INPUT NAME=\"password\" TYPE=\"hidden\" VALUE=\"%s\">\n", user->password );

	fprintf( output, "<TABLE WIDTH=\"50%%\">\n<TR>\n" );
	fprintf( output, "<TD ALIGN=\"RIGHT\">Old Password:</TD>\n" );
	fprintf( output, "<TD><INPUT TYPE=\"password\" NAME=\"oldpassword\" SIZE=\"20\"></TD>\n</TR>\n<TR>\n" );
	fprintf( output, "<TD ALIGN=\"RIGHT\">New Password:</TD>\n" );
	fprintf( output, "<TD><INPUT TYPE=\"password\" NAME=\"newpassword1\" SIZE=\"20\"></TD>\n</TR>\n<TR>\n" );
	fprintf( output, "<TD ALIGN=\"RIGHT\">Repeat New Password:</TD>\n" );
	fprintf( output, "<TD><INPUT TYPE=\"password\" NAME=\"newpassword2\" SIZE=\"20\"></TD>\n</TR>\n</TABLE>\n" );

	fprintf( output, "\n<P>\n\n" );

	fprintf( output, "<TABLE BORDER=\"0\" WIDTH=\"70%%\">\n" );
	fprintf( output, "   <TD ALIGN=\"RIGHT\" WIDTH=\"25%%\">\n" );
	fprintf( output, "   <INPUT TYPE=\"SUBMIT\" VALUE=\"Change Password\">\n" );
	fprintf( output, "   <TD ALIGN=\"RIGHT\" WIDTH=\"25%%\">\n" );
	fprintf( output, "   <INPUT TYPE=\"RESET\"  VALUE=\" Clear Fields \">\n" );
	fprintf( output, "   <TD ALIGN=\"RIGHT\" WIDTH=\"25%%\">\n" );
	fprintf( output, "   <INPUT TYPE=\"SUBMIT\" VALUE=\"   Cancel   \" NAME=\"cancel\">\n" );
	fprintf( output, "</TABLE>\n</FORM>\n\n" );

  fprintf( output, "<HR NOSHADE><P>" );

  }


void printAccountForm( FILE *output, USER *user )
  {

  fprintf( output, "<H1>Account Options</H1><P>" );
  fprintf( output, "<HR NOSHADE><P>" );

	fprintf( output, "<FORM ACTION=\"" ACCOUNT_DETAILS "\" METHOD=\"POST\">\n" );
	fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">\n", user->name );
	fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"password\" VALUE=\"%s\">\n", user->password );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"stage\" VALUE=\"1\">\n");
	fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"Change your password\">\n" );
	fprintf( output, "</FORM>\n" );
	fprintf( output, "<P>\n" );

	fprintf( output, "<FORM ACTION=\"" ACCOUNT_DETAILS "\" METHOD=\"POST\">\n" );
	fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"user\" VALUE=\"%s\">\n", user->name );
	fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"password\" VALUE=\"%s\">\n", user->password );
  fprintf( output, "<INPUT TYPE=\"hidden\" NAME=\"stage\" VALUE=\"2\">\n");
	fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"Check machine resource usage\">\n" );
	fprintf( output, "</FORM>\n" );
	fprintf( output, "<P>\n" );

	fprintf( output, "<FORM ACTION=\"" LOGIN_URL "\" METHOD=\"POST\">\n" );
	fprintf( output, "<INPUT TYPE=\"SUBMIT\" VALUE=\"Logout of system\">\n" );
	fprintf( output, "</FORM>\n" );
	fprintf( output, "<P>\n" );

  fprintf( output, "<HR NOSHADE><P>" );

  }



/* -- Program Entry Point -----------------------------------------------------
*/

int main( int argc, char *argv[] )
  {
  CGI      *cgi;
	USER     *user;
  HOSTLIST *hostList;
	char     *stagestring;
	int      stage;

  cgi = getCGIEnvironment( argc, argv );
  if( cgi == NULL )
    {
    printf( "<HR><P><B>ERROR:</B> getCGIEnvironment() failed in main()</P>"
            "<P>Email "WEBMASTER"<I>immediately</I></P>\n" );
    return 0;
    }

  user = getUser( cgi );

	hostList = newHostList();
	readHostList( hostList, ESU_MASTER_HOST_LIST );
  
  if(!(stagestring = lookupString(cgi,"stage")))
    {
    stage=0;
    }
  else
    {
    stage=atoi(stagestring);
    }


  
	borderGenerateHeader( stdout, user, False, NULL, NULL );
  
	/* Switch between the operation choices */
  switch (stage)
    {
		/* Main job control menu */
  	case 0:
			{
			printAccountForm( stdout, user );
			break;
		  }

		/* Set password */
  	case 1:
			{
			printChangeForm( stdout, user );
			break;
		  }

		/* Resources used */
  	case 2:
			{
			printResourceUsage( stdout, user, hostList );
			break;
		  }


		default:
      {
      writeError("Invalid stage number");
      break;
      }
    } 

	borderGenerateFooter( stdout );

  return 0;
  }

