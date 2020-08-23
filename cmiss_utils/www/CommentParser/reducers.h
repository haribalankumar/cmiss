/*
 *  Reducers.h
 *
 *    String reduction methods (one ofr variables, one for commands, etc)
 */

#ifndef REDUCERS_H
#define REDUCERS_H


/* -- Module Datatypes --------------------------------------------------------
*/

typedef void REDUCTION_FUNCTION( char *destination, const char *source );



/* -- Module Public Functions -------------------------------------------------
*/

void reduceVariable( char *destination, const char *source );
void reduceModule( char *destination, const char *source );
void reduceSubroutine( char *destination, const char *source );
void reduceFunction( char *destination, const char *source );
void reduceComment( char *destination, const char *source );
void reduceCommand( char *destination, const char *source );
void reduceBlockData( char *destination, const char *source );


#endif
