/*
 *  Get User.h
 *
 *    Public interface for obtaining a user from a CGI structure 
 */
#ifndef GET_USER_H
#define GET_USER_H


/* -- Required Includes -------------------------------------------------------
*/

#include "user-managment.h"
    /* For: USER */

#include "cgi-decode.h"
    /* For: CGI */



/* -- Public Method Prototypes ------------------------------------------------
*/
 
USER *getUser( CGI *cgi );
USER *getUserCrypt( CGI *cgi );
 

#endif
