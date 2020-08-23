/*
 File: cm.c
 =================
 
 This file is the main `back-end' CMISS main program.

 Copyright 1996-2006, Auckland UniServices Ltd.

*/

/* Included files */

#include <setjmp.h>
#include <signal.h>
#ifdef sun
#include <ucontext.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#if defined (unix) || defined (_AIX)
#include <sys/utsname.h>
#endif /* defined (unix) || defined (_AIX) */
#ifdef mips
#include <sys/systeminfo.h>
#  if defined(DEBUG)
#    include <libexc.h>
#  endif
/* to allow it to work on IRIX 6.2 and 6.3 */
#  define _OLD_TERMIOS 
#endif
#if defined (unix) || defined (_AIX)
#include <termios.h>
#endif /* defined (unix) || defined (_AIX) */
#include "message.h"
#include "perl_interpreter.h"

/* Defines */

#ifdef VMS
#  define CmDestroy CMDESTROY
#  define CmInitialise CMINITIALISE
#  define CmMain CMMAIN
#  define SetCmissCommandParams SETCMISSCOMMANDPARAMS
#endif
#if defined (unix) || defined (_AIX) || defined (WIN32)
/* Functions defined here referenced from Fortran */
#  define AddStringToBuffer addstringtobuffer_
#  define GetLine getline_
#  define GetSystemInfo getsysteminfo_
#  define ResetFatalHandler resetfatalhandler_
#  define SetFatalHandler setfatalhandler_
/* Fortran functions accessed */
#  define CmDestroy cmdestroy_
#  define CmInitialise cminitialise_
#  define GetNextComfileLine getnextcomfileline_
#  define GetPrevComfileLine getprevcomfileline_
#  define MainCmLoop maincmloop_
#  define PeriodicTask periodictask_
#  define SetCmissCommandParams setcmisscommandparams_
#endif


/* Defines */

#define NUMBUFFERS 20
#define LINELENGTH 255

#define NUL 0x00 /* Null character */
#define BEL 0x07 /* Bell character */
#define BS  0x08 /* Backspace character */
#define BSA 0X7f /* Alternate backspace character */
#define ESC 0x1B /* Escape character */
#define LF  0x0A /* Linefeed character */
#define SPC 0x20 /* Space character */

#define TIMEOUTTIME 1

#define MIN_PORT_NUMBER 1001
#define MAX_PORT_NUMBER 9999

/* whether we include the alignment checking code */
/* for more details, see man 3 fixade */
#define NOT_CHECK_ALIGNMENT
#if defined (CHECK_ALIGNMENT)
#  if defined (mips)
/* prototype the functions */
void handle_unaligned_traps(void);
void list_by_addr(void);
void summary_listing(void);
void print_unaligned_summary(void);
#  else /* mips */
/* not allowed unless on SGI */
#    error CHECK_ALIGNMENT not available unless using SGI
#  endif /* mips */
#endif /* CHECK_ALIGNMENT */

/* Type definitions for Fortran */

typedef int integer;

/* Function prototypes */

/* Functions defined here referenced from Fortran */
void AddStringToBuffer(char *string);
void GetLine(char *prompt, 
  char *string,
  integer *err,
  char *error_string);
void GetSystemInfo(integer size,
		   char *sysname,
		   char *winsystem);
void ResetFatalHandler(void);
void SetFatalHandler(void);

/* Fortran functions accessed */
void CmDestroy(integer *err);
void CmInitialise(integer *batch_mode,
		  integer *port1,
		  integer *port2,
		  integer *use_socket,
		  integer *err);
void GetNextComfileLine(integer *maxlength,
			char *line,
			integer *error_code);
void GetPrevComfileLine(integer *maxlength,
			char *line,
			integer *error_code);
void MainCmLoop(const integer *batch_mode,
		const char *examplenum,
		const char *comfilename,
		const char *parameterfilename,
		const integer *first,
		const char *executestring,
		integer *error_code);
void PeriodicTask(void);
void SetCmissCommandParams(const char *cmgui_arg,
			   const char *const *examplepath,
			   const char *cm_version,
			   const char *imagename,
			   const integer *idle_time,
			   integer *err);
void set_interpreter_ptr_in_comm00_(struct Interpreter **interpreter);
#ifdef _AIX /* can't find a prototype for this anywere */
void xl__sigdump(int sig,
		 int code,
		 struct sigcontext *sc);
#endif

/* Internal functions */
static void backspace(void);
static void bell(void);
static void clearline(size_t charnum,
  size_t linesize,
  char *line);
static void CmMain(integer batch_mode,
		   const char *examplenum,
		   const char *comfilename,
		   const char *parameterfilename,
		   const char *executestring,
		   integer *err);
#if defined (unix) || defined (_AIX)
static void fatalhandler(int sig,
#  if defined (sun)
                         siginfo_t *sip,
                         ucontext_t *uap);
#  else
			 int code,
			 struct sigcontext *sc);
#  endif
#endif /* defined (unix) || defined (_AIX) */
static void initfatalhandler(void);
static void putline(char *line);
static void resetterminal(void);
static void setterminal(void);
static int StringAbbrev(char *string,
  char *comstring,
  size_t minchar);
static char *strncpy0( char *dest, const char *src, size_t n );

/* Static variables */

static const char cm_version[] = "2.1";
static const char default_comfile[] = "cmiss.com";

#if defined (unix) || defined (_AIX)
static sigjmp_buf jump_buffer;
static struct sigaction fatal_sigaction;
static struct sigaction old_SIGBUS_action;
#ifdef SIGEMT
static struct sigaction old_SIGEMT_action;
#endif
/* KAT 26Feb99:
  In the debug version of CMISS the trap_uninitialized flag is on.
  This causes an SIGFPE signal every time an INVALID fpe occurs.
  Unfortunately DDOT in the sgi 7.2.1 blas library seems to perform INVALID
  operations that do not affect the result.  These are trapped.
  To override this, set the environment variable TRAP_FPE to 
  "INVALID=APPROPRIATE".  This causes the libfpe library which is linked into
  the debug version to quietly determine an appropriate result for an INVALID
  operation.
  To avoid fatalhandler interupting this process it is not used for SIGFPE.
 */
/* KAT 26Feb99:
   fatalhandler is now used for SIGFPE because many machines now have
   different blas libraries.  For the TRAP_FPE variable to have an effect
   the cmiss command 'set fatal off must be executed.
 */
static struct sigaction old_SIGFPE_action;
static struct sigaction old_SIGILL_action;
static struct sigaction old_SIGINT_action;
static struct sigaction old_SIGABRT_action;
static struct sigaction old_SIGSEGV_action;
static struct sigaction old_SIGTRAP_action;
#endif /* defined (unix) || defined (_AIX) */

static integer buffnum=0,numbuff=0;
static char buffer[NUMBUFFERS][LINELENGTH];

#if defined (unix) || defined (_AIX)
static struct termios terminal,oldterminal;
#endif /* defined (unix) || defined (_AIX) */

/* External symbols accessed */
extern const char revision_time_string[]; /* revision_time.c */

/* Code */

int main(int argc,char *argv[], char *envp[])
{
  int return_code;
  integer err=0,batch_mode=0,i,idle_time=60,
    port1=1,port2=1,show_help=0,use_socket=0;
  const char *examplepath = NULL, *cmgui_arg = "";
  const char *examplenum = "", *comfilename = "";
  const char *parameterfilename = "", *executestring = "";
  const char *imagename;
  struct stat buf;
  
  /* Initialise variables */
  return_code = EXIT_FAILURE; /* bad return code */

#if defined (mips)  
#  if defined (CHECK_ALIGNMENT)
  /* turn on alignment checking */
  handle_unaligned_traps_();
  /* list in terms of data addresses */
  list_by_addr_();
#  endif /* CHECK_ALIGNMENT */
#endif /* mips */
  
  /* Get command line arguments */

  imagename = argv[0];

  for( i=1; i < argc; i++ )
    {
      if(StringAbbrev(argv[i],"-batch",1))
	{
	  batch_mode=1;
	}
      else if(StringAbbrev(argv[i],"-cmgui",3))
      {
        if(i+1 >= argc)
	  {
	    err=1;
	    show_help=1;
	    display_message(ERROR_MESSAGE,"no cmgui argument string");
	  }
        else
	  {
	    i++;
	    cmgui_arg = argv[i];
	  }
      }
      else if(StringAbbrev(argv[i],"-epath",3))
      {
        if(i+1 < argc)
        {
          i++;
	  examplepath = argv[i];
        }
        else
        {
          err=1;
          show_help=1;
          display_message(ERROR_MESSAGE,"no example path");
        }
      }
      else if(StringAbbrev(argv[i],"-example",3))
	{
	  if(i+1 < argc)
	    {
	      i++;            
	      examplenum = argv[i];
	    }
	  else
	    {
	      err=1;
	      show_help=1;
	      display_message(ERROR_MESSAGE,"no example number");
	    }
	}
      else if(StringAbbrev(argv[i],"-execute",3))
	{
	  /* SAB 10 Oct 2000 Execute the command line string */
	  if(i+1 < argc)
	    {
	      i++;
	      executestring = argv[i];
	    }
	  else
	    {
	      err=1;
	      show_help=1;
	      printf(">>ERROR: No execution string\n");
	    }
	}
      else if(StringAbbrev(argv[i],"-h",1))
	{
	  show_help=1;
	}
      else if(StringAbbrev(argv[i],"-idle_time",3))
	{
	  if(i+1 < argc)
	    {
	      idle_time=atoi(argv[i+1]);
	      i++;
	    }
	  else
	    {
	      err=1;
	      show_help=1;
	      printf(">>ERROR: No idle time\n");
	    }
	}
      else if(StringAbbrev(argv[i],"-parameters",3))
      {
        if(i+1 < argc)
        {
          i++;            
          parameterfilename = argv[i];
        }
        else
        {
          err=1;
          show_help=1;
          display_message(ERROR_MESSAGE,"no parameters filename");
        }
      }
      else if(StringAbbrev(argv[i],"-port",3))
	{
	  if(i+1 < argc)
	    {
	      port1=atoi(argv[i+1]);
	      port2=port1+1;
	      i++;
	    }
	  else
	    {
	      err=1;
	      show_help=1;
	      printf(">>ERROR: No port number\n");
	    }
	}
      else if(StringAbbrev(argv[i],"-socket",2))
	{
	  use_socket=1;
	}
      else
	{
	  comfilename = argv[i];
	}
    }

  if(use_socket) 
  {
    if(MIN_PORT_NUMBER >= port1 || MAX_PORT_NUMBER <= port1)
    {
      err=1;
      show_help=1;
      fprintf(stderr,"main.  Invalid port number\n");
    }
    if(MIN_PORT_NUMBER >= port2 || MAX_PORT_NUMBER <= port2)
    {
      err=1;
      show_help=1;
      fprintf(stderr,"main.  Invalid port number\n");
    }
  }
  if(show_help)
  {
    printf("Usage : cm <-h>\n");
    printf("           <-batch>\n");
    printf("           <-cmgui argument_string>\n");
    printf("           <-example example_number>\n");
    printf("           <-execute execute_string>\n");
    printf("           <-epath path_to_examples_directory>\n");
    printf("           <-idle_time #>\n");
    printf("           <-socket -port #>\n");
    printf("           <-parameters parameterfilename>\n");
    printf("           <commandfilename>[cmiss]\n");
  }
  else
  {
      
    /* If no error in the command line arguments then startup CMISS */
    
    if(!err)
    {
      if(!batch_mode && getenv("CM_AUTOBATCH") && comfilename[0])
	{
	  batch_mode = 1;
	}
      if(!comfilename[0])
	{
	  if(!examplenum[0] && -1 != stat(default_comfile,&buf))
	    {
	      comfilename = default_comfile;
	    }
	}
      printf("CMISS(cm) version %s  %s\n",cm_version,revision_time_string);
      if(use_socket)
	    {
	      printf("  Using socket interface\n");
	      printf("    with ports %d and %d\n",port1,port2);
	    }
      
      /* Initialise CMISS system variables. NOTE: error trapping
      is not set here so any errors will be fatal. */
      CmInitialise(&batch_mode,&port1,&port2,&use_socket,&err);

      /* Only create interpreter if not running under cmgui */
      if(!err && !cmgui_arg[0])
      {
	 int status;
	 struct Interpreter *interpreter;
	 create_interpreter(argc,argv,comfilename,&interpreter,&status);
	 if(status)
	 {
	    interpreter_set_display_message_function(interpreter,
	       display_message, &status);
	    if (!status) err = 1;
	    set_interpreter_ptr_in_comm00_(&interpreter);
	 }
	 else
	 {
	    err = 1;
	 }
      }
      
      /* If no example path supplied then check environment */
      /* This pointer is found after create_interpreter as perl seems to move
         the environment or something? */
      if(!examplepath)
	{
	  examplepath = getenv("CMISS_EXAMPLES");
	}

      if(!err)
      {
	/* Get data from command parameters */
	SetCmissCommandParams( cmgui_arg, &examplepath,
			       cm_version, imagename, &idle_time, &err);

        if(!err)
        {
          /* Call the main back-end program*/
          CmMain(batch_mode,examplenum,comfilename,
            parameterfilename,executestring,&err);
          if(!err)
          {
            return_code = EXIT_SUCCESS; /* successful completion */
          }
          /* Call the cm destroy routine */
          CmDestroy(&err);
          if(err)
          {
            return_code = EXIT_FAILURE; /* unsuccessful completion */
            fprintf(stderr,"main.  Error calling CmDestroy\n");
          }
        }
        else
	  {
	    fprintf(stderr,
		    "main.  Error initialising command line parameters\n");
	  }
      }
      else
      {
	fprintf(stderr,"main.  Error in initialization\n");
      }
    }
    else
    {
      fprintf(stderr,"main.  Error in command line arguments\n");
    }
  }

#if defined (mips)  
#  if defined (CHECK_ALIGNMENT)
  /* print the output */
 print_unaligned_summary_();
#  endif /* CHECK_ALIGNMENT */
#endif /* mips */
 
 return return_code;
}

void AddStringToBuffer(char *string)

/*
C#### Function: AddStringToBuffer
C###  Description:
C###    Adds a string to the command buffer.
*/

{
  integer i;

  /* Place the string into the buffer if it is different from
  the last line entered */

  if(NUMBUFFERS == numbuff)
  {
    if(strcmp(buffer[NUMBUFFERS-1],string))
    {
      for(i=0; i < NUMBUFFERS-1; i++)
      {
        strcpy(buffer[i],buffer[i+1]);
      }
      strcpy(buffer[NUMBUFFERS-1],string);
    }
  }
  else
  {
    if(numbuff > 0)
    {
      if(strcmp(buffer[numbuff-1],string))
      {
        strcpy(buffer[numbuff],string);
        numbuff++;
      }
    }
    else
    {
      strcpy(buffer[numbuff],string);
      numbuff++;
    }
  }
  buffnum=numbuff;
}

void GetLine(char *prompt,
  char *string,
  integer *err,
  char *error_string)
/*
C#### Function: GetLine
C###  Description:
C###    Writes the prompt to the string and returns the line input from
C###    the terminal.
*/

{
  integer comfile_err,finished=0;
  char line[LINELENGTH],dummy[LINELENGTH];
  int intchar; /* getchar returns int */
  size_t charnum,i,promptlen,numchar;
  
  static int linelength = LINELENGTH;
    
  *err=0;
  promptlen=strlen(prompt);
  if(LINELENGTH > promptlen) 
  {
    /* Setup terminal characteristics */
    setterminal();
    /* Set the prompt */
    printf("> %s",prompt);
    strcpy(line,prompt);
    numchar=promptlen;
    charnum=promptlen;
    if(0 < promptlen)
    {
      line[promptlen]=SPC;
      numchar++;
      charnum++;
      putchar(SPC);
    }
    /* Loop until a Line Feed is found */
    while(!finished)
    {
      intchar=getchar(); /* Get the character */
      switch (intchar)
      {
        case EOF: /* Timeout has occured */
        {
          PeriodicTask();
        }; break;
        case BS: /* Backspace */
        case BSA:
        {
          if(charnum > 0) 
          {
            charnum--;
            /* Move the cursor over the backspaced character */
            putchar(BS);
            /* Adjust the characters in front of the cursor */
            for(i=charnum; i < numchar-1; i++)
            {
              line[i]=line[i+1];
              putchar(line[i]);
            }
            /* Clear the last character */
            putchar(SPC);
            /* Reposition the cursor */
            for(i=numchar; i > charnum; i--)
            {
              putchar(BS);
            }
            numchar--;
          }
          else
          {
            bell();
          }
        }; break;
        case LF: /* Newline */
        {
          finished=1;
          putchar(LF);
        }; break;
        case ESC: /* Escape character */
        {
          intchar=getchar(); /* Read [ */
          intchar=getchar(); /* Get Escape character */
          switch(intchar)
          {
            case 0x41: /* Up arrow */
            {
              if(buffnum > 0 && numbuff > 0)
              {
                /* Clear current line */
                clearline(charnum,numchar,line);
                /* Put buffer into the line */
                buffnum--;
                strcpy(line,buffer[buffnum]);
                charnum=strlen(line);
                numchar=charnum;
                /* Write out the line */
                putline(line);
              }
              else
              {
                bell();
              }
            }; break;
            case 0x42: /* Down arrow */
            {
              if(buffnum == numbuff || numbuff == 0) 
              {
                bell();
              }
              else
              {
                if(buffnum == (numbuff-1))
                {
                  /* Clear current line */
                  clearline(charnum,numchar,line);
                  /* Put prompt into the line */
                  buffnum++;
                  strcpy(line,prompt);
                  numchar=promptlen;
                  charnum=promptlen;
                  if(0 < promptlen)
                  {
                    line[promptlen]=SPC;
                    numchar++;
                    charnum++;
                  }
                  /* Write out the line */
                  putline(line);
                }
                else
                {
                  /* Clear current line */
                  clearline(charnum,numchar,line);
                  /* Put buffer into line */
                  buffnum++;
                  strcpy(line,buffer[buffnum]);
                  charnum=strlen(line);
                  numchar=charnum;
                  /* Write out the line */
                  putline(line);
                }
              }
            }; break;
            case 0x43: /* Right arrow */
            {
              if(charnum == numchar)
              {
                bell();
              }
              else
              {
                /* Move the cursor right */
                putchar(line[charnum]); 
                charnum++;
              }
            }; break;
            case 0x44: /* Right arrow */
            {
              if(charnum > 0) 
              {
                /* Move the cursor left */
                putchar(BS);
                charnum--;
              }
              else
              {
                bell();
              }
            }; break;
#ifdef mips
            case 0x31: /* extended escape sequence */
            {
              intchar=getchar(); /* Read character 0x36 */
              intchar=getchar(); /* Get identification character */
              switch (intchar)
              {
                case 0x33: /* Alt+Up arrow */
#else
            case 0x50: /* Alt+Up arrow */
#endif
                {
                  GetPrevComfileLine(&linelength,dummy,&comfile_err);
                  /* 1==got prev line ok
                  0==cant get an earlier line
                  -1==error occurred somewhere */ 
                  if(!comfile_err)
                  {
                    /* Clear current line */
                    clearline(charnum,numchar,line);
                    /* Put string into the line */
                    strcpy(line,dummy);
                    charnum=strlen(line);
                    numchar=charnum;
                    buffnum=numbuff;
                    /* Write out the line */
                    putline(line);
                  }
                  else
                  {
                    bell();
                  }
                }; break;
#ifdef mips
                case 0x36: /* Alt+Down arrow */
#else
	  case 0x51: /* Alt+Down arrow */
#endif
                {
                  GetNextComfileLine(&linelength,dummy,&comfile_err);
                  /* 1==got next line ok
                  0==no more lines to get
                  -1==error occurred somewhere */ 
                  if(!comfile_err)
                  {
                    /* Clear current line */
                    clearline(charnum,numchar,line);
                    /* Put string into the line */
                    strcpy(line,dummy);
                    charnum=strlen(line);
                    numchar=charnum;
                    buffnum=numbuff;
                    /* Write out the line */
                    putline(line);
                  }
                  else
                  {
                    bell();
                  }
                }; break;
                default:
                {
                  /* do nothing */
                }; break;
              }
#ifdef mips
              intchar=getchar(); /* Read character 0x71 */
            }; break;
            default:
            {
              /* do nothing */
            }; break;
          } /* end of escape character switch */
#endif
        }; break;
      default:
        {
          if((intchar>=SPC) && (intchar<=0x7E))
          {
            /* Insert the character */
            if(LINELENGTH-1 == numchar)
            {
              bell();
            }
            else
            {
              /* Output the character */
              putchar(intchar);
              /* Move all characters left of the cursor left */
              for(i=charnum; i < numchar; i++)
              {
                putchar(line[i]);
              }
              /* Reposition the cursor and update the line */
              for(i=numchar; i > charnum; i--)
              {
                line[i]=line[i-1];
                putchar(BS);
              }
              line[charnum]=(char)intchar;
              charnum++;
              numchar++;
            }
          }
        }; break;
	    }
    }
    /* NUL terminate the string */
    line[numchar]=NUL;
    strcpy(string,line);
    /* Reset terminal characteristics */
    resetterminal();
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: Prompt is longer than line length\n");
  }
}

void GetSystemInfo(integer size,
		   char *sysname,
		   char *winsystem)

/*
C#### Function: GetSystemInfo
C###  Description:
C###    Writes the system information of the current host
*/

{
  const char unknown[]="UNKNOWN";
  const char none[]="NONE";
#if defined (unix) || defined (_AIX)
  struct utsname name;
#endif /* defined (unix) || defined (_AIX) */
  const char *display,*nodename;

#if defined (unix) || defined (_AIX)
  if(uname(&name) < 0)
    {
      strncpy0(sysname,unknown,size);
      nodename = unknown;
    }
  else
    {
      printf(    "System nodename:   %s\n",name.nodename);
      nodename = name.nodename;
      printf(    "Machine type:      %s\n",name.machine);
#ifdef mips
      {
	const size_t numsize = 10, procsize = 200;
	char numstr[numsize], procstr[procsize];
	if(sysinfo(_MIPS_SI_NUM_PROCESSORS,numstr,numsize) <= 0)
	  {
	    strncpy0(numstr,unknown,numsize);
	  }
	if(sysinfo(_MIPS_SI_PROCESSORS,procstr,procsize) <= 0)
	  {
	    strncpy0(procstr,unknown,procsize);
	  }
	if(strcmp(numstr,"1") == 0)
	  {
	    printf("Processor type:    %s\n",procstr);
	  }
	else
	  {
	    printf("Num. processors:   %s\n",numstr);
	    printf("Processor types:   %s\n",procstr);
	  }
      }
#endif /* mips */
#ifdef _AIX /* concatenate version.release */
      printf(    "Operating system:  %s %s.%s\n",
	     name.sysname,name.version,name.release);
#else /* ! _AIX */
      printf(    "Operating system:  %s %s\n",name.sysname,name.release);
#endif /* _AIX */
      strncpy0(sysname,name.sysname,size);
    }
#else /* defined (unix) || defined (_AIX) */
  strncpy0(sysname,unknown,size);
  nodename = unknown;
#endif /* defined (unix) || defined (_AIX) */

  display=getenv("DISPLAY");
  if( !display || !*display ) 
    {
      printf(    "Display:           %s\n",none);
      strncpy0(winsystem,none,size);
    }
  else
    {      
      if(':' == display[0])
	{
	  printf("Local display:     %s%s\n",nodename,display);
	}
      else
	{
	  printf("Remote display:    %s\n",display);
	}
      strncpy0(winsystem,"MOTIF",size);
    }
}

void ResetFatalHandler()
/*
C#### Function: ResetFatalHandler
C###  Description:
C###    Resets old (default) signal handlers.
*/

{
#if defined (unix) || defined (_AIX)
#if defined (SIGBUS)
  if( 0 != sigaction(SIGBUS,&old_SIGBUS_action,NULL) )
    {
      fprintf(stderr,"WARNING: could not reset SIGBUS handler");
    }
#endif /* defined (SIGBUS) */
#ifdef SIGEMT
  if( 0 != sigaction(SIGEMT,&old_SIGEMT_action,NULL) )
    {
      fprintf(stderr,"WARNING: could not reset SIGEMT handler");
    }
#endif
  if( 0 != sigaction(SIGFPE,&old_SIGFPE_action,NULL) )
    {
      fprintf(stderr,"WARNING: could not reset SIGFPE handler");
    }
  if( 0 != sigaction(SIGILL,&old_SIGILL_action,NULL) )
    {
      fprintf(stderr,"WARNING: could not reset SIGILL handler");
    }
  if( 0 != sigaction(SIGINT,&old_SIGINT_action,NULL) )
    {
      fprintf(stderr,"WARNING: could not reset SIGINT handler");
    }
  if( 0 != sigaction(SIGABRT,&old_SIGABRT_action,NULL) )
    {
      fprintf(stderr,"WARNING: could not reset SIGABRT handler");
    }
  if( 0 != sigaction(SIGSEGV,&old_SIGSEGV_action,NULL) )
    {
      fprintf(stderr,"WARNING: could not reset SIGSEGV handler");
    }
#if defined (SIGTRAP)
  if( 0 != sigaction(SIGTRAP,&old_SIGTRAP_action,NULL) )
    {
      fprintf(stderr,"WARNING: could not reset SIGTRAP handler");
    }
#endif /* defined (SIGTRAP) */
#endif /* defined (unix) || defined (_AIX) */
}

void SetFatalHandler(void)
     
/*
C#### Function: SetFatalHandler
C###  Description:
C###    Sets the signal handlers back to be the fatalhandler.
*/

{
#if defined (unix) || defined (_AIX)
#if defined (SIGBUS)
  if( 0 != sigaction(SIGBUS,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,"WARNING: could not set SIGBUS handler");
    }
#endif /* defined (SIGBUS) */
#ifdef SIGEMT
  if( 0 != sigaction(SIGEMT,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,"WARNING: could not set SIGEMT handler");
    }
#endif
  if( 0 != sigaction(SIGFPE,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,"WARNING: could not set SIGFPE handler");
    }
  if( 0 != sigaction(SIGILL,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,"WARNING: could not set SIGILL handler");
    }
  if( 0 != sigaction(SIGINT,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,"WARNING: could not set SIGINT handler");
    }
  if( 0 != sigaction(SIGABRT,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,"WARNING: could not set SIGABRT handler");
    }
  if( 0 != sigaction(SIGSEGV,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,"WARNING: could not set SIGSEGV handler");
    }
#if defined (SIGTRAP)
  if( 0 != sigaction(SIGTRAP,&fatal_sigaction,NULL) )
    {
      fprintf(stderr,"WARNING: could not set SIGTRAP handler");
    }
#endif /* defined (SIGTRAP) */
#endif /* defined (unix) || defined (_AIX) */
}

static void CmMain(integer batch_mode,
		   const char *examplenum,
		   const char *comfilename,
		   const char *parameterfilename,
		   const char *executestring,
		   integer *err)

/*
C#### Function: CmMain
C###  Description:
C###    Calls the Main CM loop provides a place to recover from fatal signals
*/

{
  integer first=1;
  int sig_code;

#ifdef mips
  /* KAT: This environment variable affects code compiled under debug
     with SGI f77.  It makes the code send abort (ABRT) signals instead
     of just aborting when runtime errors detected in the fortran
     (e.g. index out of range).  It may also cause a core dump. */
  putenv("f77_dump_flag=y");
#endif

#if defined (unix) || defined (_AIX)
  /* Setup but don't activate fatal signal handler */
  initfatalhandler();
  /* All signals are restored (including non-fatal) if we return here
     with non-zero sig_code. */
  sig_code = sigsetjmp(jump_buffer,1);
  if(sig_code)
    {
      /* If a signal is intercepted, fatalhandler will return here. */
      /* Assuming we got here from fatalhandler, reactivate fatalhandler. */
      SetFatalHandler();

      first=0;
    }
#endif /* defined (unix) || defined (_AIX) */

  /* Enter the main CMISS command loop */
  {
    MainCmLoop( &batch_mode, examplenum, comfilename, parameterfilename,
		&first, executestring, err );
  }

}

#if defined (unix) || defined (_AIX)
static void fatalhandler(int sig,
#  if defined (sun)
                         siginfo_t *sip,
                         ucontext_t *uap)
#  else
			 int code,
			 struct sigcontext *sc)
#  endif

/*
C#### Function: fatalhandler
C###  Description:
C###    Long jumps back to CmMain for signal processing.
*/

{
#if defined(_AIX)
  /* this from libxlf90.a provides a good description of what went wrong */
  xl__sigdump(sig,code,sc);
#else
  switch(sig)
    {
#if defined (SIGBUS)
    case SIGBUS:
      {
	fprintf(stderr,">>ERROR: Bus error occured\n");
      } break;
#endif /* defined (SIGBUS) */
#if defined (SIGEMT)
    case SIGEMT:
      {
	fprintf(stderr,">>ERROR: EMT occured\n");
      } break;
#endif /* defined (SIGEMT) */
    case SIGFPE:
      {
	fprintf(stderr,">>ERROR: Floating point exception occured\n");
      } break;
    case SIGILL:
      {
	fprintf(stderr,">>ERROR: Illegal instruction occured\n");
      } break;
    case SIGINT:
      {
	fprintf(stderr,">>ERROR: Interrupt occured\n");
      } break;
    case SIGABRT:
      {
	fprintf(stderr,">>ERROR: Abort occured\n");
      } break;
    case SIGSEGV:
      {
	fprintf(stderr,">>ERROR: Segment violation occured\n");
      } break;
#if defined (SIGTRAP)
    case SIGTRAP:
      {
	fprintf(stderr,">>ERROR: Trace Trap occured\n");
      } break;
#endif /* defined (SIGTRAP) */
    default:
      {
	fprintf(stderr,">>ERROR: unrecognized signal %i occured\n", sig);
      } break;
    }

#  if defined(mips) && defined(DEBUG)
  trace_back_stack_and_print();
#  endif
#endif /* AIX */

  siglongjmp(jump_buffer,sig);
}
#endif /* defined (unix) || defined (_AIX) */

#if defined (unix) || defined (_AIX)
static void initfatalhandler(void)

/*
C#### Function: initfatalhandler
C###  Description:
C###    Initializes fatal_sigaction and records the old (default)
C###    signal handlers.
*/

{
  fatal_sigaction.sa_flags = SA_NODEFER;
  fatal_sigaction.sa_handler = (void (*)(int))fatalhandler;
  if( 0 != sigemptyset(&fatal_sigaction.sa_mask) )
    {
      fprintf(stderr,"WARNING: sigemptyset failed in initfatalhandler");
    }

#if defined (SIGBUS)
  sigaction(SIGBUS,NULL,&old_SIGBUS_action);
#endif /* defined (SIGBUS) */
#if defined (SIGEMT)
  sigaction(SIGEMT,NULL,&old_SIGEMT_action);
#endif /* defined (SIGEMT) */
  sigaction(SIGFPE,NULL,&old_SIGFPE_action);
  sigaction(SIGILL,NULL,&old_SIGILL_action);
  sigaction(SIGINT,NULL,&old_SIGINT_action);
  sigaction(SIGABRT,NULL,&old_SIGABRT_action);
  sigaction(SIGSEGV,NULL,&old_SIGSEGV_action);
#if defined (SIGTRAP)
  sigaction(SIGTRAP,NULL,&old_SIGTRAP_action);
#endif /* defined (SIGTRAP) */
}
#endif /* defined (unix) || defined (_AIX) */

static int StringAbbrev(char *string,
  char *comstring,
  size_t minchar)

/*
C#### Function: StringAbbrev
C###  Description:
C###    Returns 1 if the string is an abbreviation (with a minimum
C###    match of minchar) of comstring, 0 if not.
*/

{
  if(strlen(string) <= strlen(comstring))
  {
    if(strlen(string) >= minchar)
    {
      if(!strncmp(string,comstring,strlen(string)))
      {
        return 1;
      }
      else
      {
        return 0;
      }
    }
    else
    {
      return 0;
    }
  }
  else
  {
    return 0;
  }
}

static void backspace(void)
{
  putchar(BS);
  putchar(SPC);
  putchar(BS);
}

static void bell(void)
{
  putchar(BEL);
}

static void clearline(size_t charnum,
  size_t linesize,
  char *line)
{
  size_t i;

  for(i=charnum; i < linesize; i++)
  {
    putchar(SPC);
  }
  for(i=linesize; i > 0; i--)
  {
    backspace();
    line[i]=NUL;
  }
  line[0]=NUL;
}

static void putline(char *line)
{
  size_t linesize,i;

  linesize=strlen(line);
  for(i=0;i < linesize; i++)
  {
    putchar(line[i]);
  }
}

static void resetterminal(void)
{
#if defined (unix) || defined (_AIX)
  /* cpb use POSIX termios instead of termio */
  tcsetattr(0, TCSANOW, &oldterminal);
#endif /* defined (unix) || defined (_AIX) */
}

static void setterminal(void)
{
#if defined (unix) || defined (_AIX)
  /* cpb use POSIX termios instead of termio */
  tcgetattr(0, &terminal);
  tcgetattr(0, &oldterminal);
  /* Set terminal echo off */
  terminal.c_lflag &= ~ECHO;
  /* Set mode to non-canocial input */
  terminal.c_lflag &= ~ICANON;
  /* Set non-canocial input parameters */
  terminal.c_cc[VMIN] = 0;
  terminal.c_cc[VTIME] = TIMEOUTTIME;
  /* cpb use POSIX termios instead of termio */
  tcsetattr(0, TCSANOW, &terminal);
#endif /* defined (unix) || defined (_AIX) */
}

static char *strncpy0( char *dest, const char *src, size_t n )

/*
  Perform a function similar to strncpy but ensure that the
  destination dest is null terminated.  n is the number of bytes
  available in dest not the number of characters to copy.  i.e. it
  includes the null terminator.
*/

{
  *dest = NUL;
  return strncat( dest, src, n - 1 );
}

/*
    Local Variables: 
    tab-width: 8
    End: 
*/
