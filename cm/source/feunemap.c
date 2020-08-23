/*
  Wrapper routines for cm to interface with unemap routines
  Created : Leo Cheng
  Date    : June 2002
*/

#include <stdarg.h>
#include <string.h>
#include "debug.h"
#include "feunemap.h"
#include "message.h"

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
  int end_search,start_search;

  int return_code=0;
  int number_of_events=1,objective_values_step=1;
  float gain=1.0, offset=0.0;

  int *event_times_addr,i;


  return_code=calculate_time_series_objective
    (
      *detection_type,
      *objective_type,
      *width,
      gain,
      offset,
      *num_samples,
      objective_values_step,
      objective_values);
  
  /*check return codes if == 0 then error occurred*/
  start_search=0;
  end_search=*num_samples-1;

  if(return_code == 0)
  {
    *err=1;
    strcpy(error_string,
      ">>ERROR: return_code error calculate_time_series_objective");   
  }
  else
  {
    return_code=calculate_time_series_event_markers(
      start_search,
      end_search,            
      *detection_type,
      objective_values,
      *num_samples,
      objective_values_step,
      number_of_events,
      *threshold_percentage,
      *minimum_separation_milliseconds,
      *level,
      *frequency,
      events_found,      
      &event_times_addr
      );

    if(return_code == 0)
    {
      *err=1;
      strcpy(error_string,
        ">>ERROR: return_code error in calculate_time_series_event_markers");   
    }
    else
    {
      
      /*Copy the answers back into the static arrays*/
      if(*events_found > *max_events)
      {
        printf("WARNING: Found %3d events. More than can be stored, increase NEVENTSM in code\n",*events_found);
        *events_found=*max_events;
      }
      for(i=0;i<*events_found;i++)
      {
        *(events+i)=*(event_times_addr+i);
      }
    }
  
  }  
}

