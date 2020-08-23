/*
  Wrapper routines for cm to interface with unemap routines
  when the libunemap.a library cannot be found
  Created : Leo Cheng
  Date    : June 2002
*/

#include <stdarg.h>
#include "feunemap.h"

/********************************/

void CmUnemapWrapper(
  enum Event_detection_algorithm *detection_type,
  int *events,
  int *events_found,
  int *max_events,
  int *minimum_separation_milliseconds,
  int *num_samples,
  int *num_signals,
  enum Event_detection_objective *objective_type,
  int *width,
  int *threshold_percentage,
  float *objective_values,
  float *level,
  float *frequency,
  int *err,
  char *error_string)
{

  *err=1;
  strcpy(error_string,
    ">>ERROR: Need to link with unemap library to use event detection");   
  
}

