/*
 *  CMISS Line IO.h
 *
 *    CMISS derived class of the Line IO base class
 */
#ifndef CMISS_LINEIO_H
#define CMISS_LINEIO_H


/* -- Required Include Files --------------------------------------------------
*/

#include "lineio.h"
    /* For: FILESTATE, newVirtualFileState() */



/* -- Constructors ------------------------------------------------------------
*/

FILESTATE *newCmissFileState( char *fileName );


#endif
