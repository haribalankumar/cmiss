/*
 *  Rule Parse.h
 * 
 *    The Rule Parsing modules public interface
 */

#ifndef RULE_PARSE_H
#define RULE_PARSE_H


/* -- Required Includes -------------------------------------------------------
*/

#include "rules.h"
    /* For: RULES */



/* -- Public Methods ----------------------------------------------------------
*/

int parseRuleFile( RULES *rules, const char *ruleFileName );


#endif
