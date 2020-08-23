/*******************************************************************************
FILE : write-files.c

LAST MODIFIED : 6 January 1997

DESCRIPTION :
This program reads files from the standard input, and writes them out.
==============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* #include "<project includes specified relative to project base directory eg. */
/*         graphics/finite_element.h> */

/*
Global variables
----------------
*/
/* <variable declarations.  Same as in .h except without extern and with */
/*         initializers> */

/*
Module types
------------
*/
enum Input_state
/*******************************************************************************
LAST MODIFIED : 6 January 1997

DESCRIPTION :
This type is used to determine which state the program is in.
NOFILE - Waiting to recieve a file.
INFILE - Recieving and writing a file.
==============================================================================*/
{
  STATE_NOFILE,
  STATE_INFILE
}; /* enum Input_state */

/*
Module variables
----------------
*/
/* <declarations for variables which are global to this module only.  Should be */
/*         storage class static> */

/*
Module functions 
----------------
*/
/* <definitions for functions which are used in this module only.  Should be */
/*         storage class static> */

/*
Global functions
----------------
*/
/* <definitions for the functions prototyped in .h> */
#define MAX_LINE_LENGTH 1000

void main()
{
  char filename[MAX_LINE_LENGTH],line[MAX_LINE_LENGTH];
  FILE *file;
  enum Input_state state;

  
  state=STATE_NOFILE;
  while(!feof(stdin))
  {
    line[MAX_LINE_LENGTH-1] = '\0';
    if((1==fscanf(stdin,"%[^\n]",line))&&(0==fscanf(stdin,"\n")))
    {
      /* check for the flag to make sure we haven't overwritten */
      if(!line[MAX_LINE_LENGTH-1])
      {
        switch(state)
        {
          case STATE_NOFILE:
          {
            if(!strncmp("START:",line,6))
            {
              strcpy(filename,&(line[6]));
/*             strcat(filename,"new");*/
              if((file=fopen(filename,"w"))!=NULL)
              {
                state=STATE_INFILE;
              }
              else
              {
                fprintf(stderr,"Could not open the file %s\n",filename);
              }
            }
            else
            {
              fprintf(stderr,"Incorrect tag\n");
            }
          } break;
          case STATE_INFILE:
          {
            if(!strncmp("DATA:",line,5))
            {
              strcpy(line,&(line[5]));
              fprintf(file,"%s\n",line);
            }
            else if (!strncmp("END:",line,4))
            {
              if(!fclose(file))
              {
                state=STATE_NOFILE;
              }
              else
              {
                fprintf(stderr,"Could not close the file\n");
              }
            }
            else 
            {
              fprintf(stderr,"Incorrect tag\n");
            }
          } break;
          default:
          {
            fprintf(stderr,"Invalid state\n");
          } break;
        }
      }
      else
      {
        fprintf(stderr,"line too long - recompile\n");
      }
    }
    else
    {
      fprintf(stderr,"Could not read in line\n");
    }
  }
}
