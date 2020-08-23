/*
 *  Password.h
 *
 *    Public interface for editing passwords 
 */
#ifndef PASSWORD_H
#define PASSWORD_H


/* -- Required Includes -------------------------------------------------------
*/


/* -- Public Method Prototypes ------------------------------------------------
*/
 
int changePassword( char *user, char *password );
int addUser( char *user, char *password, char *email, char *fullName, char *home );

#endif
