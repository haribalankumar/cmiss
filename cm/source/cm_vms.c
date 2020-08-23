/*
File: cm_vms.c
=================
 
This file provides c system dependendant routines for VMS.

Functions included:
 
(externally referenced)
CmMain             Main cm routine (calls main loop)

*/

/* Included files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Defines */

#define CmMain CMMAIN
#define MainCmLoop MAINCMLOOP

/* Type definitions */

typedef int integer;

/* Function prototypes */

/* External functions */
void CmMain(integer example,
  char *examplenum,
  char *comfilname,
  integer parameters,
  char *parameterfilename);
void MainCmLoop(integer *example,
  char *examplenum,
  char *comfilename,
  integer *parameters,
  char *parameterfilename,
  integer *first);
/* Internal functions */

/* Global variables */

/* Code */

void CmMain(integer example,
  char *examplenum,
  char *comfilename,
  integer parameters,
  char *parameterfilename)

/*
C#### Function: CmMain
C###  Description:
C###    Calls the Main CM loop
*/

{
  integer first=0;
  
  /* Enter the main CMISS command loop */
  
  MainCmLoop(&example,examplenum,comfilename,&parameters,
    parameterfilename,&first);

}





