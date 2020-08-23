/*******************************************************************************
FILE : example_path.c

LAST MODIFIED : 8 October 2002

DESCRIPTION :
==============================================================================*/
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <signal.h>
#include <sys/wait.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include "debug.h"
#include "message.h"
#include "example_path.h"

/*
Global functions
----------------
*/

char *resolve_example_path(char *example_path, char *directory_name,
	char **comfile_name, char **requirements)
/*******************************************************************************
LAST MODIFIED : 17 April 2000

DESCRIPTION :
Uses the executable $example_path/common/resolve_example_path to demangle
a short example name into a full path.  The returned string is ALLOCATED.
Code is basically a repeat from general/child_process but don't want to create
an object and want to encapsulate the whole process in one file for calling
from the back end.
<*comfile_name> is allocated and returned as well if the resolve function
returns a string for it.  This too must be DEALLOCATED by the calling function.
<*requirements> is either set to NULL or is an allocated string specifying 
the features required to run this example, i.e. whether the example 
needs cmgui and/or cm. The requirements are comma separated.  This too must 
be DEALLOCATED by the calling function.
==============================================================================*/
{
	char *return_string = NULL;
#if defined (unix) || defined (_AIX)
#define BLOCKSIZE (100)
	char *filename, last_char, *new_string,
	  *comfile_offset, *end_offset, *requirements_offset;
	fd_set readfds;
	int status, stdin_filedes, stdout_filedes[2];
	pid_t process_id;
	size_t index, string_size;
	ssize_t number_read;
	struct stat buf;
	struct timeval timeout_struct;
#endif /* defined (unix) || defined (_AIX) */

	ENTER(resolve_example_path);

	if (example_path && directory_name)
	{
#if defined (unix) || defined (_AIX)
		if (ALLOCATE(filename, char, strlen(example_path) + 
			strlen(directory_name) +50))
		{
			sprintf(filename, "%s/common/resolve_example_path", example_path);

			if (-1 == stat(filename,&buf))
			  {
				display_message( ERROR_MESSAGE,
								 "File %s does not exist", filename);
			  }
 
			else if ( -1 == pipe(stdout_filedes) )
			  {
				display_message (ERROR_MESSAGE,
"resolve_example_path. Unable to create pipe: %s", strerror(errno) );
			  }
			else
			  {
			    process_id = fork();

				if (0 == process_id)
				  {
					/* Child process comes here */

					close(stdout_filedes[0]);
				
					/* The child shouldn't read anything */
					/* Is this the best way to redirect stdin to /dev/null? */
					stdin_filedes = open ("/dev/null", O_RDONLY);
					dup2 (stdin_filedes,STDIN_FILENO);
					close(stdin_filedes);
					/* Remap stdout */
					dup2 (stdout_filedes[1],STDOUT_FILENO);
					close(stdout_filedes[1]);

					/* Execute the filename */
/* !!! Should first ensure that all non-system file descriptors are closed! */
					execlp (filename, filename, directory_name, (char *)0);
					/* execlp only returns on error
					   (as on success the process gets overlayed). */
					display_message (ERROR_MESSAGE, "%s: %s",
									 filename, strerror(errno) );
					exit(EXIT_FAILURE);
				  }

				/* Parent (or no fork) */
				close(stdout_filedes[1]); /* This was for the child. */

				if( -1 == process_id )
				  {
					display_message (ERROR_MESSAGE,
"resolve_example_path.  Unable to fork process: %s", strerror(errno) );
					close(stdout_filedes[0]);
				  }
				else /* Have child */
				  {
					/* We get a SIGPIPE if we write when the process has
					   close stdin or exited.  Either handle the SIGPIPE or,
					   since resolve_example_path doesn't read from stdin,
					   don't write.
					*/
					/* 					    write(stdin_filedes[1], directory_name, strlen(directory_name) + 1); */

					FD_ZERO(&readfds);
					FD_SET(stdout_filedes[0], &readfds);
					string_size = 2 * BLOCKSIZE;
					last_char = 0xff;
					index = 0;
					if (return_string = ALLOCATE(return_string, char,
												 2 * BLOCKSIZE))
					  {
						while((last_char != 0) && return_string)
						  {
							if (index + BLOCKSIZE > string_size)
							  {
								if(REALLOCATE(new_string, return_string, char,
											  index + 2 * BLOCKSIZE))
								  {
									return_string = new_string;
									string_size = index + 2 * BLOCKSIZE;
								  }
								else
								  {
									display_message(ERROR_MESSAGE,"resolve_example_path."
													"  Unable to reallocate string");
									DEALLOCATE(return_string);
									return_string = (char *)NULL;
								  }
							  }
							if (return_string)
							  {
								int select_code;
								timeout_struct.tv_sec = 20;
								timeout_struct.tv_usec = 0;
								/* !!! Select returns -1 for errors including EINTR */
								do
								  {
									select_code =
									  select( FD_SETSIZE, &readfds,
											  NULL, NULL,	&timeout_struct );
								  }
								while( select_code == -1 && errno == EINTR );

								if( select_code > 0 )
								  {
									while( (number_read =
											read( stdout_filedes[0],
												  (void *)
												  (return_string + index),
												  BLOCKSIZE )
											) == -1
										   && errno == EINTR )
									  ;
									if (/* error */
										-1 == number_read
									    /* end of file before terminating 0 */
										|| 0 == number_read
										)
									  {
										display_message (ERROR_MESSAGE,
"Error reading from resolve_example_path: %s",
number_read ? strerror(errno) : "Incomplete" );
										DEALLOCATE(return_string);
										return_string = (char *)NULL;
									  }
									else
									  /* We may not read the \n from the
										 child.  Hopefully it doesn't matter
										 if it gets a SIGPIPE. */
									  while((last_char != 0) && number_read)
										{
										  last_char = *(return_string + index);
										  number_read--;
										  index++;
										}
								  }
								else
								  {
									if( select_code == 0 )
									  {
										display_message(ERROR_MESSAGE,
"resolve_example_path.  Timed out waiting for response from child");
									  }
									else
									  {
										display_message
										  ( ERROR_MESSAGE, "select failed: %s", strerror(errno) );
									  }
									DEALLOCATE(return_string);
									return_string = (char *)NULL;
								  }
							  }
						  }
					  }
					else
					  {
						display_message(ERROR_MESSAGE,"resolve_example_path."
										"  Unable to allocate string");
					  }
/* 						write(stdin_filedes[1], end, 3); */
/* 					    close(stdin_filedes[1]); */

					close(stdout_filedes[0]);

					if (return_string)
						{
						   comfile_offset = (char *)NULL;
						   requirements_offset = (char *)NULL;
						   /* Look for the first space separator in the 
							  returned string */
						   if (comfile_offset = strchr(return_string, ' '))
							{
							  /* Terminate the example path string */
							  *comfile_offset = 0;
							  comfile_offset++;

							  /* Look for the next space */
							  if (requirements_offset = strchr(comfile_offset, ' '))
							  {
								 /* Terminate the comfile string */
								 *requirements_offset = 0;
								 requirements_offset++;

								 /* Look for the end of this word */
								 if (end_offset = strchr(requirements_offset, ' '))
								 {
									/* Terminate the requirements string */
									*end_offset = 0;
								 }
							  }
						   }
						   if (comfile_name)
						   {
							  if (comfile_offset && ALLOCATE(*comfile_name, char,
								 strlen (comfile_offset) + 1))
							  {
								 strcpy (*comfile_name, comfile_offset);
							}
							else
							{
								*comfile_name = (char *)NULL;
							}
						   }
						   if (requirements)
						   {
							  if (requirements_offset && ALLOCATE(*requirements, 
								 char, strlen (requirements_offset) + 1))
							  {
								 strcpy (*requirements, requirements_offset);
							  }
							  else
							  {
								 *requirements = (char *)NULL;
							  }
						   }
							if (ALLOCATE(new_string, char,
								strlen(return_string) + strlen(example_path) + 5))
							{
								sprintf(new_string, "%s/%s/", example_path,
									return_string);
								DEALLOCATE(return_string);
								return_string = new_string;
							}
							else
							{
								DEALLOCATE(return_string);
								display_message(ERROR_MESSAGE,"resolve_example_path."
									"  Unable to make final reallocate of string");
								return_string = (char *)NULL;
							}
						}

					/* Reap the child */

					/* Assuming that if we got what we expected from the
					   child, then it is probably going to exit.  Is this
					   reasonable? */
					if ( ! return_string )
					  {
						/* Otherwise tell the child to finish.  If it won't
						   accept a SIGTERM then I assume it has a good reason
						   not to exit yet? */
						kill (process_id, SIGTERM);
					  }
					waitpid (process_id, &status, 0);
				  }
			  }

		}
		else
		{
			display_message(ERROR_MESSAGE,"resolve_example_path. Unable to allocate program string");
		}
#else /* defined (unix) || defined (_AIX) */
		display_message(ERROR_MESSAGE,"resolve_example_path.  "
			"Not implemented yet.");
#endif /* defined (unix) || defined (_AIX) */
	}
	else
	{
		display_message(ERROR_MESSAGE,"resolve_example_path.  Invalid argument(s).");
	}

	LEAVE;

	return (return_string);
} /* resolve_example_path */

void destroy_example_path(char **example_path, char**comfile_name)
/*******************************************************************************
CREATED : KAT 19 April 2000
LAST MODIFIED : 19 April 2000

DESCRIPTION :
Destroys the example path created by resolve_example_path.
==============================================================================*/
{
  ENTER(destroy_example_path);

  if(example_path && comfile_name)
    {
      if(*example_path)
	{
	  DEALLOCATE(*example_path);
	  *example_path = NULL;
	}
      else
	{
	  display_message(WARNING_MESSAGE,"destroy_example_path.  Example path does not exist");
	}
      if(*comfile_name)
	{
	  DEALLOCATE(*comfile_name);
	  *comfile_name = NULL;
	}
    }
  else
    {
      display_message(ERROR_MESSAGE,"destroy_example_path.  Invalid argument(s).");
    }

  LEAVE;
}

/* Local Variables:  */
/* tab-width: 4 */
/* End:  */
