/*
File: fecutils.c
=================
 
This file provides c utility routines for cmiss.

Functions included:

BinaryCloseFile      Close a binary file
BinaryOpenFile       Open a binary file
BinaryReadFile       Read data from a binary file
BinarySetFile        Sets the position of a binary file
BinarySkipFile       Skip bytes in a binary file
BinaryWriteFile      Write data to a binary file
CFreeTimer           Free a timer for reuse
CGetTimer            Allocate a new timer
CResetTimers         Resets (frees) all timer handles
C_getenv             Returns a pointer to a c string containing the value of
                     an environment variable
C_strlen             Returns the length of a c string
CStringLen           Sets a variable to the length of a c string
CTimer_CPU           Performs CPU timing
CTimer_REAL          Performs real timing
FreeMemory           Deallocates dynamic memory
Get_Seconds          Returns the number of seconds since a fixed time
IsBinaryFileOpen     Returns whether or not a binary file is open
IsEndBinaryFile      Returns whether or not at eof of a binary file
MallocMemory         Allocates dynamic memory
PackCharacters       Packs characters into a c string
SystemCommand        Issues a command to the operating system
UnpackCharacters     Unpacks characters from a c string

*/

/* Included files */

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef VMS
#include <time.h>
#include <types.h>

#elif POSIX_TIMERS
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>

#elif BSD_TIMERS
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

#else /* WIN32/ANSI_C_TIMERS */
#include <time.h>

#endif

#include "message.h"

/* Defines */

#ifdef VMS
#define BinaryCloseFile BINARYCLOSEFILE
#define BinaryOpenFile BINARYOPENFILE
#define BinaryReadFile BINARYREADFILE
#define BinarySetFile BINARYSETFILE
#define BinarySkipFile BINARYSKIPFILE
#define BinaryWriteFile BINARYWRITEFILE
#define CFreeTimer CFREETIMER
#define CGetTimer CGETTIMER
#define CResetTimers CRESETTIMERS
#define C_getenv C_GETENV
#define C_strlen C_STRLEN
#define CStringLen CSTRINGLEN
#define CTimer_CPU CTIMER_CPU
#define CTimer_REAL CTIMER_REAL
#define FreeMemory FREEMEMORY
#define Get_Seconds GET_SECONDS
#define IsBinaryFileOpen ISBINARYFILEOPEN
#define IsEndBinaryFile ISENDBINARYFILE
#define MallocMemory MALLOCMEMORY
#define BigMallocMemory BIGMALLOCMEMORY
#define PackCharacters PACKCHARACTERS
#define SystemCommand SYSTEMCOMMAND
#define UnPackCharacters UNPACKCHARACTERS
#endif
#if defined(unix) || defined (_AIX) || defined (WIN32)
#define BinaryCloseFile binaryclosefile_
#define BinaryOpenFile binaryopenfile_
#define BinaryReadFile binaryreadfile_
#define BinarySetFile binarysetfile_
#define BinarySkipFile binaryskipfile_
#define BinaryWriteFile binarywritefile_
#define CFreeTimer cfreetimer_
#define CGetTimer cgettimer_
#define CResetTimers cresettimers_
#define C_getenv c_getenv_
#define C_strlen c_strlen_
#define CStringLen cstringlen_
#define CTimer_CPU ctimer_cpu_
#define CTimer_REAL ctimer_real_
#define FreeMemory freememory_
#define Get_Seconds get_seconds_
#define IsBinaryFileOpen isbinaryfileopen_
#define IsEndBinaryFile isendbinaryfile_
#define MallocMemory mallocmemory_
#define BigMallocMemory bigmallocmemory_
#define ReallocMemory reallocmemory_
#define PackCharacters packcharacters_
#define Unlink_File unlink_file_
#define SystemCommand systemcommand_
#define UnPackCharacters unpackcharacters_
#endif


/* Item Type defines */
/*   These should be the same as those in mach00.cmn */
#define INTEGERTYPE  1
#define FLOATTYPE    2
#define DOUBLETYPE   3
#define CHARTYPE     4
#define LOGICALTYPE  5
#define SHORTINTTYPE 6
#define COMPLEXTYPE  7
#define DOUBLECOMPLEXTYPE 8 
#define POINTERTYPE 9 

/* Binary file defines */
#define MAXBINFILES 99

/* Type definitions */

typedef int integer;
typedef integer logical;
struct Time_info
/* This structure holds the zero values for cpu and elapsed time */
{
  long cpu,elapsed;    /* zero values */
  integer in_use;        /* is timer currently in use */
};

/* Function prototypes */

void BinaryCloseFile(integer *fileid,
  integer *err, 
  char *error_string);
void BinaryOpenFile(integer *fileid,
  char *filename,
  char* access_code,
  integer *err,
  char *error_string);
void BinaryReadFile(integer *fileid,
  integer *endian,
  integer *number_of_items, 
  integer *item_type,
  char *data,
  integer *err,
  char *error_string);
void BinarySetFile(integer *fileid,
  integer *set_code,
  integer *err,
  char *error_string);
void BinarySkipFile(integer *fileid,
  integer *number_of_bytes, 
  integer *err,
  char *error_string);
void BinaryWriteFile(integer *fileid,
  integer *endian,
  integer *number_of_items, 
  integer *item_type,
  char *data,
  integer *err,
  char *error_string);
char *C_getenv(const char *name);
integer C_strlen(const char *string);
void CStringLen(int *length,
  char *string);
void CFreeTimer(integer *ihandle,
  integer *err,
  char *error_string);
void CGetTimer(integer *ihandle,
  integer *err,
  char *error_string);
void CResetTimers(integer *err,
  char *error_string);
void CTimer_CPU(double *return_time,
  integer *flag,
  integer *err,
  char *error_string);
void CTimer_REAL(double *return_time,
  integer *flag,
  integer *err,
  char *error_string);
void FreeMemory(void **ptr,
  integer *err,
  char *error_string);
void Get_Seconds(integer *num_seconds_ptr);
void IsBinaryFileOpen(int *fileid,
  int *returncode,
  int *err,
  char *error_string);
void IsEndBinaryFile(int *fileid,
  int *returncode,
  int *err,
  char *error_string);
void MallocMemory(integer *number_of_items,
  integer *item_type,
  integer *init_flag,
  void **ptr,
  integer *err,
  char *error_string,
  integer error_size );
void BigMallocMemory(integer *number_of_items1,
  integer *number_of_items2,
  integer *item_type,
  integer *init_flag,
  void **ptr,
  integer *err,
  char *error_string,
  integer error_size );
void ReallocMemory(integer *number_of_items,
  integer *item_type,
  void **ptr,
  integer *err,
  char *error_string,
  integer error_size );
void PackCharacters(integer *integer_char,
  integer *char_num,
  char *integer_string);
void Unlink_File( const char *path, integer *err );
void SystemCommand(char *command,
  integer *err,
  char *error_string); 
void UnPackCharacters(integer *integer_char,
  integer *char_num,
  char *integer_string);

void iload(integer *n, integer *const, integer *x, integer *incx);
void sload(integer *n, double *const, double *x, integer *incx);

/* Global variables */

FILE *binaryfiles[MAXBINFILES];
#define MAX_TIMERS 10
/* WARNING!!! If you change MAX_TIMERS, make sure that there are enough
initialisers for timer_information!!! */
static struct Time_info time_information[MAX_TIMERS] =
{
  {0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},
  {0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}
};

/* Code */


static char *strncpy0( char *dest, const char *src, size_t n )

/*
  Perform a function similar to strncpy but ensure that the
  destination dest is null terminated.  n is the number of bytes
  available in dest not the number of characters to copy.  i.e. it
  includes the null terminator.
*/

{
  *dest = 0;
  return strncat( dest, src, n - 1 );
}

void BinaryCloseFile(integer *fileid,
  integer *err,
  char *error_string)

/*
C#### Function: BinaryCloseFile
C###  Description:
C###    Closes the binary file specified by fileid.
*/

{
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    if(!binaryfiles[*fileid-1])
    {
      *err = 0;
      /* cpb 31/5/95 Don't signal an error when closing a file that
      is not open 
      *err=1;
      strcpy(error_string,">>ERROR: binary file is not open"); */
    }
    else
    {
      *err = fclose(binaryfiles[*fileid-1]);
      binaryfiles[*fileid-1] = (FILE *)NULL;
      if(*err != 0)
      {
        strcpy(error_string,">>ERROR: error closing binary file");
      }
    }
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}

void BinaryOpenFile(integer *fileid,
  char *filename,
  char *access_code,
  integer *err,
  char *error_string)

/* 
C#### Function: BinaryOpenFile
C###  Description:
C###    Opens a binary file specified by fileid and name filename.
*/

{
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    if(binaryfiles[*fileid-1])
    {
      *err=1;
      strcpy(error_string,">>ERROR: binary file is already open");
    }
    else
    {
      if(binaryfiles[*fileid-1] = fopen(filename,access_code))
	    {
	      *err=0;
	    }
      else
	    {
	      *err=1;
	      strcpy(error_string,">>ERROR: binary file could not be opened");
	    }
    }
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}

void BinaryReadFile(integer *fileid,
  integer *endian,
  integer *number_of_items, 
  integer *item_type,
  char *data,
  integer *err,
  char *error_string)

/* 
C#### Function: BinaryReadFile
C###  Description:
C###    Reads number_of_items of data of a type given by item_type from
C###    a binary file specified by fileid into an array iven by data.

The default endian ordering is big endian. This is specified by
endian=0. If little endian is required endian must be set to a number
other than 0.

*/

{
  integer i,item_size,j,number_of_bytes,start_byte,temp;
  FILE* binaryfile;
  
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    binaryfile=binaryfiles[*fileid-1];
    if(!binaryfile)
    {
      *err = 1;
      strcpy(error_string,">>ERROR: binary file is not open");
    }
    else
    {
      if(INTEGERTYPE == *item_type)
	    {
	      /* Integer data */
	      item_size=sizeof(integer);
	    }
      else if(FLOATTYPE == *item_type)
	    {
	      /* Float data */
	      item_size=sizeof(float);
	    }
      else if(DOUBLETYPE == *item_type)
	    /* Double data */
	    {
	      item_size=sizeof(double);
	    }
      else if(CHARTYPE == *item_type)
	    /* Character data */
	    {
	      item_size=sizeof(char);
	    }
      else if(LOGICALTYPE == *item_type)
	    /* Logical data */
	    {
	      item_size=sizeof(logical);
	    }
      else if(SHORTINTTYPE == *item_type)
	    /* Short integer data */
	    {
	      item_size=sizeof(short int);
	    }
      else
	    {
	      *err=1;
	      strcpy(error_string,">>ERROR: Invalid item type");
	    }
	  
      if(0 == *endian || CHARTYPE == *item_type) 
	    {
	      /* Default big endian format */
	      number_of_bytes=*number_of_items * item_size;
	      for(i = 0; i < number_of_bytes ; i++)
        {
          temp = getc(binaryfile); 
          data[i] = (char)temp;  
        }
        if(CHARTYPE == *item_type)
        {
          data[number_of_bytes]='\0';
        }
	    }
      else
	    {
	      /* Little endian format - must reverse byte ordering */
	      for(i = 0; i < *number_of_items; i++)
        {
          start_byte=i*item_size;
          for(j = 0; j < item_size; j++)
          {
            temp = getc(binaryfile); 
            data[start_byte+item_size-j-1] = (char)temp;
          }
        }
	    }
	  
      *err=ferror(binaryfile);
      if(*err != 0)
	    {
	      strcpy(error_string,">>ERROR: error reading binary file");
	    }
    }
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}

void BinarySetFile(integer *fileid,
  integer *set_code,
  integer *err,
  char *error_string)

/* 
C#### Function: BinarySetFile
C###  Description:
C###    Sets the position of the file pointer of the binary file (given
C###    by fileid) to either the beginning (set_code=0), current
C####   position (set_code=1) or end (set_code=2) of the file.
*/

{
  FILE* binaryfile;
  
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    binaryfile=binaryfiles[*fileid-1];
    if(!binaryfile)
    {
      *err = 1;
      strcpy(error_string,">>ERROR: binary file is not open");
    }
    else
    {
      switch(*set_code)
      {
        case 0: /* Beginning of a file */
        {
          *err=fseek(binaryfile,(long)0,SEEK_SET);
          if(*err != 0)
          {
	    display_message( WARNING_MESSAGE, "Could not set to beginning of file: %s",
			     strerror(errno) );
            strcpy(error_string,">>ERROR: Could not set to beginning of file");
          }
        } break;
        case 1: /* Current file position */
        {
          *err=fseek(binaryfile,(long)0,SEEK_CUR);
          if(*err != 0)
          {
	    display_message( WARNING_MESSAGE, "Could not set to current file position: %s",
			     strerror(errno) );
            strcpy(error_string,">>ERROR: Could not set to current file position");
          }
        } break;
        case 2: /* End of file */
        {
          *err=fseek(binaryfile,(long)0,SEEK_END);
          if(*err != 0)
          {
	    display_message( WARNING_MESSAGE, "Could not set to end of file: %s",
			     strerror(errno) );
            strcpy(error_string,">>ERROR: Could not set to end of file position");
          }
        } break;
        default:
        {
          *err=1;
          strcpy(error_string,">>ERROR: Invalid set_code");          
        } break;
      }
    }
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}

void BinarySkipFile(integer *fileid,
  integer *number_of_bytes, 
  integer *err,
  char *error_string)

/* 
C#### Function: BinarySkipFile
C###  Description:
C###    Skips number_of_bytes bytes of data in a binary file specified
C###    by fileid.
*/

{
  integer i;
  FILE* binaryfile;
  
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    binaryfile=binaryfiles[*fileid-1];
    if(!binaryfile)
    {
      *err = 1;
      strcpy(error_string,">>ERROR: binary file is not open");
    }
    else
    {
      for(i = 0; i < *number_of_bytes ; i++)
	{
	  getc(binaryfile); 
	}
	  
      *err=ferror(binaryfile);
      if(*err != 0)
	    {
	      strcpy(error_string,">>ERROR: error skipping binary file");
	    }
    }
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}

void BinaryWriteFile(integer *fileid,
  integer *endian,
  integer *number_of_items, 
  integer *item_type,
  char *data,
  integer *err,
  char *error_string)

/* 
C#### Function: BinaryWriteFile
C###  Description: 
C###    Writes number_of_items of data of a type given byitem_type to a
C###    binary file specified by fileid from an array given by data. 

The default endian ordering is big endian. This is specified by
endian=0. If little endian is required endian must be set to a number
other than 0.

*/

{
  integer i,item_size,j,number_of_bytes,start_byte;
  FILE* binaryfile;
  
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    binaryfile=binaryfiles[*fileid-1];
    if(!binaryfile)
    {
      *err = 1;
      strcpy(error_string,">>ERROR: binary file is not open");
    }
    else
    {
      if(INTEGERTYPE == *item_type)
	    {
	      /* Integer data */
	      item_size=sizeof(integer);
	    }
      else if(FLOATTYPE == *item_type)
	    {
	      /* Float data */
	      item_size=sizeof(float);
	    }
      else if(DOUBLETYPE == *item_type)
	    /* Double data */
	    {
	      item_size=sizeof(double);
	    }
      else if(CHARTYPE == *item_type)
	    /* Character data */
	    {
	      item_size=sizeof(char);
	    }
      else if(LOGICALTYPE == *item_type)
	    /* Logical data */
	    {
	      item_size=sizeof(logical);
	    }
      else if(SHORTINTTYPE == *item_type)
	    /* Short integer data */
	    {
	      item_size=sizeof(short int);
	    }
      else
	    {
	      *err=1;
	      strcpy(error_string,">>ERROR: Invalid item type");
	    }
	  
      if(0 == *endian || CHARTYPE == *item_type) 
	    {
	      /* Default big endian format */
	      number_of_bytes=*number_of_items * item_size;
	      for(i = 0; i < number_of_bytes ; i++)
        {
          putc(data[i],binaryfile);
        }
	    }
      else
	    {
	      /* Little endian format - must reverse byte ordering */
	      for(i = 0; i < *number_of_items; i++)
        {
          start_byte=i*item_size;
          for(j = 0; j < item_size; j++)
          {
            putc(data[start_byte+item_size-j-1],binaryfile);
          }
        }
	    }
      
      *err=ferror(binaryfile);
      if(*err != 0)
	    {
	      strcpy(error_string,">>ERROR: error writing binary file");
	    }
    }
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}


void CFreeTimer(integer *ihandle,
  integer *err,
  char *error_string)
/*
C#### Function: CFreeTimer
C###  Description:
C###    Frees up a timer handle that was created with GetTimer
*/
{
  /* checking is performed in FREETIMER on ihandle, check again here */
  if((*ihandle>=0)&&(*ihandle<MAX_TIMERS))
  {
    if(time_information[*ihandle].in_use)
    {
      *err=0;
      time_information[*ihandle].in_use = 0;
    }
    else
    {
      *err=1;
      strcpy(error_string,">>ERROR: Timer handle not used");
    }
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: Invalid timer handle");
  }
  *ihandle = -1; /* ensure it cannot be used later */
}


void CGetTimer(integer *ihandle,
  integer *err,
  char *error_string)
/*
C#### Function: CGetTimer
C###  Description:
C###    Gets a new timer handle for use with Timer.  The handle
C###    must be freed after use via FreeTimer
*/
{
  int found,i;
  
  /* checking is performed in GETTIMER on ihandle, check again here */
  found = 0;
  for(i=0;i<MAX_TIMERS;i++)
  {
    if (!time_information[i].in_use)
    {
      found = 1;
      *ihandle = i;
    }
  }
  if(found)
  {
    *err=0;
    time_information[*ihandle].in_use = 1;
  }
  else
  {
    *ihandle = -1; /* get some attention! */
    *err=1;
    strcpy(error_string,
      ">>ERROR: No available timer handles (increase MAX_TIMERS)");
  }
}


void CResetTimers(integer *err,
  char *error_string)
/*
C#### Function: CResetTimers
C###  Description:
C###    Resets (frees) all timer handles for later use
*/
{
  int i;
  
  /* checking is performed in CFreeTimer on ihandle */
  for(i=0;i<MAX_TIMERS;i++)
  {
    if (time_information[i].in_use)
    {
      CFreeTimer(&i,err,error_string);
    }
  }
}

char *C_getenv(const char *name)
{
  return(getenv(name));
}

integer C_strlen(const char *string)
{
  return(strlen(string));
}

void CStringLen(integer *length,
  char *string)
{
  *length=strlen(string);
}


void CTimer_CPU(double *return_time,
  integer *flag,
  integer *err,
  char *error_string)
/*
C#### Function: CTimer_CPU
C###  Description:
C###    IF FLAG=CPU_USER the CPU time used while executing
C###    instructions in the user space of the calling process
C###    is returned in seconds since initialisation.
C###    IF FLAG=CPU_SYSTEM the CPU time used by the system
C###    on behalf of the calling process is returned in seconds
C###    since initialisation.
C###    IF FLAG=CPU_TOTAL the sum of the CPU time used while
C###    executing instructions in the user space of the
C###    calling process and CPU time used by the system
C###    on behalf of the calling process is returned in seconds
C###    since initialisation.
C###    IF FLAG=CPU_TICKS the routine returns the resolution of
C###    the cpu clock in ticks per second.
C###    Called by CPU_TIMER in fe01.f.
*/
/* keep these up to date with time02.cmn!!! */
#define CPU_TOTAL    1
#define CPU_USER     2
#define CPU_SYSTEM   3
#define CPU_TICKS    4
{
  double system_time,user_time,ticks;
#ifdef VMS
  struct tbuffer current_time;
#elif POSIX_TIMERS
  struct tms current_time;
#elif BSD_TIMERS
  struct rusage r;
#else /* WIN32, ANSI_C_TIMERS */
  clock_t current_time;
#endif

#ifdef VMS
  times(&current_time);

#elif POSIX_TIMERS
  times(&current_time);

#elif BSD_TIMERS
#if defined(RUSAGE_THREAD) && defined(_OPENMP)
  (void) getrusage(RUSAGE_THREAD, &r);
#else
  (void) getrusage(RUSAGE_SELF, &r);
#endif

#else /* WIN32, ANSI_C_TIMERS */
  current_time = clock();
  if (current_time == -1) {
    current_time = 0;
  }
#endif

#ifdef VMS
  ticks       = 100.0;
  user_time   = (double)current_time.proc_user_time/ticks;
  system_time = (double)current_time.proc_system_time/ticks;

#elif POSIX_TIMERS
  ticks       = sysconf(_SC_CLK_TCK);
  user_time   = (double)current_time.tms_utime/ticks;
  system_time = (double)current_time.tms_stime/ticks;

#elif BSD_TIMERS
  ticks       = 1000000.0;
  user_time   = r.ru_utime.tv_sec + (double)r.ru_utime.tv_usec/ticks;
  system_time = r.ru_stime.tv_sec + (double)r.ru_stime.tv_usec/ticks;

#else /* WIN32, ANSI_C_TIMERS */
  ticks       = CLOCKS_PER_SEC;
  user_time   = ((double)current_time)/ticks;
  system_time = 0.0;
#endif

  *err=0;
  switch(*flag)
  {
    case CPU_TOTAL:
    {
      *return_time = user_time + system_time;
    }; break;
    case CPU_USER:
    {
      *return_time = user_time;
    }; break;
    case CPU_SYSTEM:
    {
      *return_time = system_time;
    }; break;
    case CPU_TICKS:
    {
      *return_time = ticks;
    }; break;
    default:
    {
      *err=1;
      strcpy(error_string,">>ERROR: Invalid operation code");
      *return_time = -99999.0; /* get some attention! */
    }
  }     
} /* CTimer_CPU */


void CTimer_REAL(double *return_time,
  integer *flag,
  integer *err,
  char *error_string)
/*
C#### Function: CTimer_REAL
C###  Description:
C###     Returns the value of time in seconds since
C###     00:00:00 UTC, January 1,1970.
C###     Called by REAL_TIMER in fe01.f.
*/
/* keep these up to date with time02.cmn!!! */
#define REAL_TOTAL    1
{
#ifdef BSD_TIMERS
  struct timeval  tv;
  struct timezone tz;
  double ticks;
  double elapsed_time;
#else /* VMS, POSIX_TIMERS, WIN32, ANSI_C_TIMERS */
	time_t elapsed_time;
#endif

#ifdef BSD_TIMERS
  (void) gettimeofday(&tv, &tz);

  ticks        = 1000000.0;
  elapsed_time = ((double)tv.tv_sec) + ((double)tv.tv_usec)/ticks;

#else /* VMS, POSIX_TIMERS, WIN32, ANSI_C_TIMERS */
  elapsed_time = time(NULL);
#endif

  *err=0;
  switch(*flag)
  {
    case REAL_TOTAL:
    {
      *return_time = (double)elapsed_time;
    }; break;
    default:
    {
      *err=1;
      strcpy(error_string,">>ERROR: Invalid operation code");
      *return_time = -99999.0; /* get some attention! */
    }
  }     
} /* CTimer_REAL */


void FreeMemory(void **ptr,
  integer *err,
  char *error_string)

/*
C#### Function: FreeMemory
C###  Description:
C###    Deallocates dynamically allocated memory.
*/

{
  if(*ptr)
  {
    *err=0;
    free(*ptr);
    *ptr=NULL;
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: pointer to free is NULL");
  }
}

void Get_Seconds(integer *num_seconds_ptr)

/*
C#### Function: FreeMemory
C###  Description:
C###    Deallocates dynamically allocated memory.
*/

{
  *num_seconds_ptr = time(NULL);
}

void IsBinaryFileOpen(integer *fileid,
  integer *returncode,
  integer *err,
  char *error_string)

/* 
C#### Function: IsBinaryFileOpen
C###  Description:
C###    Returns 1 in returncode if a binary file specified  by fileid
C###    is open, 0 if not.
*/

{
  FILE* binaryfile;
  
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    binaryfile=binaryfiles[*fileid-1];
    if(!binaryfile)
    {      
      *returncode = 0;
    }
    else
    {
      *returncode = 1;
    }
    *err=0;
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}


void IsEndBinaryFile(integer *fileid,
  integer *returncode,
  integer *err,
  char *error_string)

/* 
C#### Function: IsEndBinaryFile
C###  Description:
C###    Returns 1 in returncode if a binary file specified by fileid
C###    is at end of file (eof), 0 if not.
*/

{
  FILE* binaryfile;
  
  if((*fileid >= 1) && (*fileid <= MAXBINFILES-1))
  {
    binaryfile=binaryfiles[*fileid-1];
    if(feof(binaryfile))
    {      
      *returncode = 1;
    }
    else
    {
      *returncode = 0;
    }
    *err=0;
  }
  else
  {
    *err=1;
    strcpy(error_string,">>ERROR: file ID is out of range");
  }
}


void MallocMemory(integer *number_of_items,
  integer *item_type,
  integer *init_flag,
  void **ptr,
  integer *err,
  char *error_string,
  integer error_size)

/*
C#### Function: MallocMemory
C###  Description:
C###    Allocates dynamic memory.
*/

{
  size_t item_size;
  
  if(*ptr)
  {
    *ptr=(void *)NULL;
    *err=1;
    strncpy0(error_string,"pointer to allocate is not NULL",
	     error_size);
  }
  else
  {
    if(0 != *number_of_items)
    {
      item_size=0;
      switch(*item_type)
      {
        case INTEGERTYPE:
        {
          /* Integer data */
          item_size=sizeof(integer);
        } break;
        case FLOATTYPE:
        {
          /* Float data */
          item_size=sizeof(float);
        } break;
        case DOUBLETYPE:
        {
          /* Double data */
          item_size=sizeof(double);
        } break;
        case CHARTYPE:
        {
          /* Character data */
          item_size=sizeof(char);
        } break;
        case LOGICALTYPE:
        {
          /* Logical data */
          item_size=sizeof(logical);
        } break;
        case SHORTINTTYPE:
        {
          /* Short integer data */
          item_size=sizeof(short);
        } break;
        case COMPLEXTYPE:
        {
          /* Complex data */
          item_size=2*sizeof(float);
        } break;
        case DOUBLECOMPLEXTYPE:
        {
          /* Double precision complex data */
          item_size=2*sizeof(double);
        } break;
        default:
        {
          item_size=0;
          *err=1;
          strncpy0(error_string,"Invalid item type",error_size);
        } break;
      }
      if(0 == *err)
      {
	if(*number_of_items < 0)
          {
            *err=1;
            strncpy0(error_string,"number_of_bytes < 0, Integer overflow?",
		     error_size);
          }
	else
	  {
	    if(1 == *init_flag)
	      *ptr = calloc(*number_of_items,item_size);
	    else
	      *ptr = malloc((size_t)*number_of_items * item_size);

	    if( ! *ptr )
	      {
		int length;
		*err=1;
		strncpy0(error_string,"Memory allocation failed: ",error_size);

		length = strlen(error_string);
		strncpy0( error_string + length, strerror(errno),
			  error_size - 1 - length );
	      }
	  }
      }
    }
    else
    {
      *ptr=(void *)NULL;
    }
  }
}

void BigMallocMemory(integer *number_of_items1,
  integer *number_of_items2,
  integer *item_type,
  integer *init_flag,
  void **ptr,
  integer *err,
  char *error_string,
  integer error_size)

/*
C#### Function: BigMallocMemory
C###  Description:
C###    Allocates dynamic memory. The number of items is split into two parts
C###    to account for possible 32 bit overflow on some platforms when 
C###    requesting a large amount of memory in a 64 bit executable. Note 
C###    that there may still be some limitations if malloc is used (as
C###    opposed to calloc) as the product of number_of_items and item_size
C###    is placed into one size_t variable.
C###  Written by: Mark Trew 17 August 2005
*/

{
  size_t item_size;
  long int number_of_items;

  number_of_items = (long int)(*number_of_items1) * (long int)(*number_of_items2);

  if(*ptr)
  {
    *ptr=(void *)NULL;
    *err=1;
    strncpy0(error_string,"pointer to allocate is not NULL",
	     error_size);
  }
  else
  {
    if(0 != number_of_items)
    {
      item_size=0;
      switch(*item_type)
      {
        case INTEGERTYPE:
        {
          /* Integer data */
          item_size=sizeof(integer);
        } break;
        case FLOATTYPE:
        {
          /* Float data */
          item_size=sizeof(float);
        } break;
        case DOUBLETYPE:
        {
          /* Double data */
          item_size=sizeof(double);
        } break;
        case CHARTYPE:
        {
          /* Character data */
          item_size=sizeof(char);
        } break;
        case LOGICALTYPE:
        {
          /* Logical data */
          item_size=sizeof(logical);
        } break;
        case SHORTINTTYPE:
        {
          /* Short integer data */
          item_size=sizeof(short);
        } break;
        case COMPLEXTYPE:
        {
          /* Complex data */
          item_size=2*sizeof(float);
        } break;
        case DOUBLECOMPLEXTYPE:
        {
          /* Double precision complex data */
          item_size=2*sizeof(double);
        } break;
        default:
        {
          item_size=0;
          *err=1;
          strncpy0(error_string,"Invalid item type",error_size);
        } break;
      }
      if(0 == *err)
      {
	if(number_of_items < 0)
          {
            *err=1;
            strncpy0(error_string,"number_of_bytes < 0, Integer overflow?",
		     error_size);
          }
	else
	  {
	    if(1 == *init_flag)
	      *ptr = calloc((size_t)number_of_items,item_size);
	    else
	      *ptr = malloc((size_t)number_of_items * item_size);

	    if( ! *ptr )
	      {
		int length;
		*err=1;
		strncpy0(error_string,"Memory allocation failed: ",error_size);

		length = strlen(error_string);
		strncpy0( error_string + length, strerror(errno),
			  error_size - 1 - length );
	      }
	  }
      }
    }
    else
    {
      *ptr=(void *)NULL;
    }
  }
}

void ReallocMemory(integer *number_of_items,
  integer *item_type,
  void **ptr,
  integer *err,
  char *error_string,
  integer error_size)

/*
C#### Function: ReallocMemory
C###  Description:
C###    Reallocates dynamic memory.
*/

{
  size_t item_size;
  void *new_ptr;
  
      item_size=0;
      switch(*item_type)
      {
        case INTEGERTYPE:
        {
          /* Integer data */
          item_size=sizeof(integer);
        } break;
        case FLOATTYPE:
        {
          /* Float data */
          item_size=sizeof(float);
        } break;
        case DOUBLETYPE:
        {
          /* Double data */
          item_size=sizeof(double);
        } break;
        case CHARTYPE:
        {
          /* Character data */
          item_size=sizeof(char);
        } break;
        case LOGICALTYPE:
        {
          /* Logical data */
          item_size=sizeof(logical);
        } break;
        case SHORTINTTYPE:
        {
          /* Short integer data */
          item_size=sizeof(short);
        } break;
        case COMPLEXTYPE:
        {
          /* Complex data */
          item_size=2*sizeof(float);
        } break;
        case DOUBLECOMPLEXTYPE:
        {
          /* Double precision complex data */
          item_size=2*sizeof(double);
        } break;
        case POINTERTYPE:
        {
          /* Pointer data */
          item_size=sizeof(void*);
        } break;
        default:
        {
          item_size=0;
          *err=1;
          strncpy0(error_string,"Invalid item type",error_size);
        } break;
      }
      if(0 == *err)
      {
	      if(*number_of_items < 0)
          {
            *err=1;
            strncpy0(error_string,"number_of_bytes < 0, Integer overflow?",
		          error_size);
          }
	      else
	      {
	        new_ptr = realloc(*ptr,(size_t)*number_of_items*item_size);

	       if( new_ptr )
         {
           *ptr=new_ptr;
         }
         else
         
         
         
         
	      {
		int length;
		*err=1;
		strncpy0(error_string,"Memory allocation failed: ",error_size);

		length = strlen(error_string);
		strncpy0( error_string + length, strerror(errno),
			  error_size - 1 - length );
	      }
	  }
      }
}

void PackCharacters(integer *integer_char,
  integer *char_num,
  char *integer_string)

/*
C#### Function: PackCharacters
C###  Description:
C###    Packs a fortran string (as a string of integers) into a c
C###    string (as a string of characters).
*/

{
  char chr;

  chr = (char)*integer_char;
  *(integer_string+*char_num) = chr;

}

void Unlink_File( const char *path, integer *err )

/*
C#### Function: Unlink_File
C###  Description: Unlinks a file.  If successful *err is set to 0.
C###    Otherwise, a warning message is produced and err is set to errno.
*/

{
  *err = 0;
  if( 0 != unlink(path) )
    {
      *err = errno;
      display_message( WARNING_MESSAGE, "Can't unlink %s: %s",
		       path, strerror(errno) );
    }
}
    
void SystemCommand(char *command,
  integer *err,
  char *error_string)

/*
C#### Function: SystemCommand
C###  Description: 
C###    Issues a command to the operating system
*/

{
  integer returncode;

  returncode = system(command);
  if(returncode)
  {
    *err=1;
    strcpy(error_string,">>ERROR: ");
    strcat(error_string,strerror(errno));
    errno=0;
  }
  else
  {
    *err=0;
  }
}


void UnPackCharacters(integer *integer_char,
  integer *char_num,
  char *integer_string)

/*
C#### Function: UnpackCharacters
C###  Description: 
C###    Unpacks a c string (as a string of characters)
C###    into a fortran string (as a string of integers).
*/

{
  char chr;

  chr = *(integer_string+*char_num);
  *integer_char = (int)chr;

}
