/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 */
#include <sys/types.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/resource.h>

double extract_slu(tv)
struct timeval *tv;
{
  double tmp;

  tmp = tv->tv_sec;
  tmp += tv->tv_usec/1000000.0;
 
  return(tmp);
}

double dclock()
{
    struct timeval tp;
/*     struct timezone tzp; */

    /* wall-clock time */
    gettimeofday(&tp,(struct timezone *)0);

    return(extract_slu(&tp));
}

