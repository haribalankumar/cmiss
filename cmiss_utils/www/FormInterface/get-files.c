/*******************************************************************************
FILE : get-files.c

LAST MODIFIED : 11 January 2002

DESCRIPTION :
This program outputs files specified using the standard input.
==============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
/* #include "<project includes specified relative to project base directory eg. " */
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
/* <declarations for types which are only used within the module.  Should be */
/*         storage class static> */

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
#define NUM_BYTES 5
#define O_RDONLY 0
#define O_WRONLY 1

void main()
{
  int fildes,readResult,filedesout;
  char filename[MAX_LINE_LENGTH];
  void *buf[MAX_LINE_LENGTH];
  FILE *file;
  size_t nbyte;
  mode_t mode = 666;

  nbyte=NUM_BYTES;
  while(!feof(stdin))
  {
    filename[MAX_LINE_LENGTH-1] = '\0';    
    if((1==fscanf(stdin,"%[^\n]",filename))&&(0==fscanf(stdin,"\n"))) 
    {
      /* check for the flag to make sure we haven't overwritten */
      if(!filename[MAX_LINE_LENGTH-1])
      {
        /* open the file */
        if((fildes=open(filename,O_RDONLY))!=-1)
        {
          if((file=fdopen(fildes,"r"))!=NULL)
          {
            fprintf(stdout,"START:%s\n",filename);
            /* Create an output file */
            if((filedesout=creat(
              "/product/cmiss/cmiss_utils/www/FormInterface/outputfile"),mode)==-1)
            {
              perror("create"); 
              exit(1); 
            }
            outfile=fdopen(filedesout,"w");

            readResult=1;
            while (readResult>0)
            {              
              /* read in  NUM_BYTES */
              readResult=read(fildes, buf, nbyte );
              
              if ( readResult > 0)
              {
                fprintf(stderr,"read result : %d\n",readResult);
                if (write(filedesout, buf, nbyte) < 0) 
                { 
                  perror("write"); 
                  exit(1); 
                } 
/*                 } */
/*                 else */
/*                 { */
/*                   fprintf(stderr,"line too long - recompile\n"); */
/*                 } */
              }    
              /*             else */
              /*             { */
              /*               fprintf(stderr,"Could not read in line\n"); */
              /*             } */
              else if( readResult < 0)
              {
                perror("read");
                exit(1);
              }
              
            }
            fprintf(stdout,"END:\n");
            if(fclose(file))
            {
              perror("fclose");
              exit(1);
            }
          }
          else
          {
            perror("fdopen");
            exit(1);
          }
        }
        else
        {
          perror("open");
          exit(1);
        }
      }
      else
      {
        fprintf(stderr,"Filename too long - recompile\n");
      }
    }
    else
    {
      fprintf(stderr,"Could not read in filename\n");
    }
  }
}


