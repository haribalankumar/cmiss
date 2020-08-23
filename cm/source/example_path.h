/*******************************************************************************
FILE : example_path.h

LAST MODIFIED : 17 April 2000

DESCRIPTION :
==============================================================================*/
#if !defined (EXAMPLE_PATH_H)
#define EXAMPLE_PATH_H

#define resolve_example_path resolve_example_path_
#define destroy_example_path destroy_example_path_

/*
Global functions
----------------
*/

char *resolve_example_path(char *example_path, char *directory_name,
	char **comfile_name, char **requirements);
/*******************************************************************************
LAST MODIFIED : 17 April 2000

DESCRIPTION :
Uses the executable $example_path/common/resolve_example_path to demangle
a short example name into a full path.  The returned string is ALLOCATED.
<*comfile_name> is allocated and returned as well if the resolve function
returns a string for it.  This too must be DEALLOCATED by the calling function.
<*requirements> is either set to NULL or is an allocated string specifying 
the features required to run this example, i.e. whether the example 
needs cmgui and/or cm. The requirements are comma separated.  This too must 
be DEALLOCATED by the calling function.
==============================================================================*/
void destroy_example_path(char **example_path, char **comfile_name);
/*******************************************************************************
CREATED : KAT 19 April 2000
LAST MODIFIED : 19 April 2000

DESCRIPTION :
Destroys the example path created by resolve_example_path.
==============================================================================*/

#endif /* !defined (EXAMPLE_PATH_H) */
