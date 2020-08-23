/* The document root for the web server. */
#define CMISS_WWW_ROOT "/hpc/cmissweb/htdocs"
/* URLPATH indicates the component of the path that must be specified in the
   url. */
#define CMISS_WWW_URLPATH "/development"
#define CM_URLPATH CMISS_WWW_URLPATH "/help/cm"
#define LOOKUP_URLPATH CM_URLPATH "/lookup"
#define COMMANDS_URLPATH CM_URLPATH "/commands"
#define EXAMPLES_URLPATH CMISS_WWW_URLPATH "/examples"

#define LOOKUP_PROGRAM "cmiss-lookup.cgi"
#define LOOKUP_PROGRAM_URLPATH "/cgi-bin/" LOOKUP_PROGRAM

