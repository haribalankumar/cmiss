/*
 *  Test Line IO.h
 *
 *    Test derived class of the Line IO base class
 */
#ifndef TEST_LINEIO_H
#define TEST_LINEIO_H


/* -- Required Include Files --------------------------------------------------
*/

#include "lineio.h"
    /* For: FILESTATE, newVirtualFileState() */



/* -- Constructors ------------------------------------------------------------
*/

FILESTATE *newTestFileState( char *fileName );


#endif
