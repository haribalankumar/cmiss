/*******************************************************************************
FILE : fzwormole.c

LAST MODIFIED : 2 October 1997

DESCRIPTION : 

Set of stub routines for inclusion in cm when the wormhole library is
not available or included.

==============================================================================*/

/* Included files */

#include <stdio.h>

#if defined (unix) || defined (_AIX) || defined (WIN32)
/* name mappings */
#define Wh_output_f_create wh_output_f_create_
#define Wh_output_f_destroy wh_output_f_destroy_
#define Wh_output_f_can_open wh_output_f_can_open_
#define Wh_output_f_can_close wh_output_f_can_close_
#define Wh_output_f_how_many_items wh_output_f_how_many_items_
#define Wh_output_f_num_items wh_output_f_num_items_
#define Wh_output_f_open_message wh_output_f_open_message_
#define Wh_output_f_close_message wh_output_f_close_message_
#define Wh_output_f_get_int wh_output_f_get_int_
#define Wh_output_f_get_double wh_output_f_get_double_
#define Wh_output_f_get_char wh_output_f_get_char_
#define Wh_output_f_get_remainder wh_output_f_get_remainder_
#define Wh_output_f_update wh_output_f_update_
#define Wh_input_f_create wh_input_f_create_
#define Wh_input_f_destroy wh_input_f_destroy_
#define Wh_input_f_open_message wh_input_f_open_message_
#define Wh_input_f_close_message wh_input_f_close_message_
#define Wh_input_f_add_int wh_input_f_add_int_
#define Wh_input_f_add_double wh_input_f_add_double_
#define Wh_input_f_add_char wh_input_f_add_char_
#define Wh_input_f_update wh_input_f_update_
#endif /* defined (unix) */


/*
output
------
*/
void Wh_output_f_create()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_output_f_create");
}

void Wh_output_f_destroy()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_output_f_destroy");
}

void Wh_output_f_can_open()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_output_f_can_open");
}

void Wh_output_f_can_close()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_output_f_can_close");
}

void Wh_output_f_num_items()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_output_f_num_items");
}

void Wh_output_f_open_message()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_output_f_open_message");
}

void Wh_output_f_close_message()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_output_f_close_message");
}


void Wh_output_f_get_int()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_output_f_get_int");
}


void Wh_output_f_get_double()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_output_f_get_double");
}


void Wh_output_f_get_char()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_output_f_get_char");
}


void Wh_output_f_get_remainder()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_output_f_get_remainder");
}


void Wh_output_f_update()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_output_f_update");
}


/*
input
------
*/

void Wh_input_f_create()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_input_f_create");
}

void Wh_input_f_destroy()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_input_f_destroy");
}

void Wh_input_f_open_message()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_input_f_open_message");
}

void Wh_input_f_close_message()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_input_f_close_message");
}

void Wh_input_f_add_int()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_input_f_add_int");
}

void Wh_input_f_add_double()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_input_f_add_double");
}

void Wh_input_f_add_char()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_input_f_add_char");
}

void Wh_input_f_update()
{
  fprintf(stderr,"ERROR: link with wormhole library: need Wh_input_f_update");
}
