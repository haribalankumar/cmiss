
/*----------------------------------------------------------------------------
 *  timer.c -- timers for the soltest code.
 *
 *  Written By: S.E. Norris
 *
 *  CPP options:
 *  BSD_TIMERS   -- use high res BSD getrusage() and gettimeofday() routines
 *  POSIX_TIMERS -- use low res POSIX times() and time() routines
 *  (default)    -- use low res ANSI C clock() and time() routines
 *----------------------------------------------------------------------------*/

#ifdef POSIX_TIMERS
#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#elif BSD_TIMERS
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#else /* WIN32, ANSI_C_TIMERS */
#include <time.h>
#endif

#include <stdio.h>

#define Timer_start     timer_start_
#define Timer_stop      timer_stop_
#define Timer_clear     timer_clear_
#define Timer_read_busy timer_read_busy_
#define Timer_read_idle timer_read_idle_
#define Timer_read_wall timer_read_wall_

#define NCLOCKS 8

struct clockdat {
  int state;
  int et[2];
  int et_old[2];
  int ct[2];
  int ct_old[2];
  int it[2];
} s[NCLOCKS] = {
  {0,{0,0},{0,0},{0,0},{0,0},{0,0}},
  {0,{0,0},{0,0},{0,0},{0,0},{0,0}},
  {0,{0,0},{0,0},{0,0},{0,0},{0,0}},
  {0,{0,0},{0,0},{0,0},{0,0},{0,0}},
  {0,{0,0},{0,0},{0,0},{0,0},{0,0}},
  {0,{0,0},{0,0},{0,0},{0,0},{0,0}},
  {0,{0,0},{0,0},{0,0},{0,0},{0,0}},
  {0,{0,0},{0,0},{0,0},{0,0},{0,0}}
};


/*----------------------------------------------------------------------------
 * Private routines
 *----------------------------------------------------------------------------*/

/*
 * CPU and wall clocks
 */
void cpu_clock(int *secs, int *usecs)
{
#ifdef POSIX_TIMERS
  /* POSIX timers */
  struct tms r;
  double t;
  double ticks;

  ticks = (double)( (clock_t) sysconf(_SC_CLK_TCK) );
  if (ticks == -1.0) {
    fprintf(stderr,"cpu_clock: error in getting CLK_TCK\n");
    exit(-1);
  }

  (void) times(&r);
  t = (double)( r.tms_utime + r.tms_stime );

  *secs  = (int) ( t/ticks );
  *usecs = (int) ( t - ( (double)(*secs)*ticks ) )*( 1.0e6/ticks );

#elif BSD_TIMERS
  /* BSD timers */
  struct rusage r;

#ifdef RUSAGE_THREAD
  /* For AIX */
  (void) getrusage(RUSAGE_THREAD, &r);
#else
  /* For IRIX, Linux, OSF, BSD */
  (void) getrusage(RUSAGE_SELF, &r);
#endif

  *secs  = (int) (r.ru_utime.tv_sec  + r.ru_stime.tv_sec); 
  *usecs = (int) (r.ru_utime.tv_usec + r.ru_stime.tv_usec);

#else /* WIN32, ANSI_C_TIMERS */
  double t;
  double ticks = (double)CLOCKS_PER_SEC;

  t = (double)clock();
  if (t == -1.0) {
    t = 0.0;
  }

  *secs  = (int) ( t/ticks );
  *usecs = (int) ( t - ( (double)(*secs)*ticks ) )*( 1.0e6/ticks );
#endif
}

void wall_clock(int *secs,int *usecs)
{
#ifdef POSIX_TIMERS
  /* POSIX timers */
  time_t t;

  t = time(NULL);

  *secs  = (int) t;
  *usecs = 0;

#elif BSD_TIMERS
  /* BSD timers */
  struct timeval  tv;
  struct timezone tz;

  (void) gettimeofday(&tv, &tz);

  *secs  = (int) tv.tv_sec; 
  *usecs = (int) tv.tv_usec;

#else /* WIN32, ANSI_C_TIMERS */
  /* ANSI C timers */
  time_t t;

  t = time(NULL);

  *secs  = (int) t;
  *usecs = 0;

#endif
}


/*
 *  Update clock values
 */
static void update_clock(int et[2], int et_old[2], int et_new[2])
{
  int tmp;

  et[0] += et_new[0] - et_old[0];
  tmp = et_new[1] - et_old[1];
  while (tmp < 0) {
    tmp += 1000000;
    et[0] -= 1;
  }
  et[1] += tmp;
  while(et[1] > 1000000) {
    et[1] -= 1000000;
    et[0] += 1;
  }

  et_old[0] = et_new[0];
  et_old[1] = et_new[1];
}

static void update_idle(int it[2], int et[2], int ct[2])
{
  int tmp;

  it[0] = et[0] - ct[0];
  tmp = et[1] - ct[1];
  while (tmp < 0) {
    tmp += 1000000;
    it[0] -= 1;
  }
  it[1] = tmp;
  while(it[1] > 1000000) {
    it[1] -= 1000000;
    it[0] += 1;
  }
}



/*---------------------------------------------------------
 * Public interface
 *---------------------------------------------------------*/

/*
 *  Start the clock
 */
int Timer_start(int *clock)
{
  int i = *clock - 1;

  if (i >= NCLOCKS || i < 0) {
    return 1;
  }

  s[i].state = 1;

  /* Start the CPU time counters */
  cpu_clock(&s[i].ct_old[0], &s[i].ct_old[1]);

  /* Start the elapsed time counters */
  wall_clock(&s[i].et_old[0], &s[i].et_old[1]);

  return 0;
}


/*
 *  Stop the clock
 */
int Timer_stop(int *clock)
{
  int i = *clock - 1;
  int ct_new[2], et_new[2];

  if (i >= NCLOCKS || i < 0) {
    return 1;
  }

  if (s[i].state == 0) {
    return 1;
  }

  s[i].state = 0;

  /* Stop the CPU counters */
  cpu_clock(&ct_new[0], &ct_new[1]);

  /* Stop the elapsed time counters */
  wall_clock(&et_new[0], &et_new[1]);

  /* Update the clock values */
  update_clock(s[i].et, s[i].et_old, et_new);
  update_clock(s[i].ct, s[i].ct_old, ct_new);
  update_idle(s[i].it, s[i].et, s[i].ct);

  return 0;
}


int Timer_clear(int *clock)
{
  int i = *clock - 1;

  if (i >= NCLOCKS || i < 0) {
    return 1;
  }

  s[i].state = 0;
  s[i].et[0] = 0;
  s[i].et[1] = 0;
  s[i].ct[0] = 0;
  s[i].ct[1] = 0;
  s[i].it[0] = 0;
  s[i].it[1] = 0;

  return 0;
}

double Timer_read_wall(int *clock)
{
  int i = *clock - 1;

  if (i < 0 || i >= NCLOCKS) {
    return 0.0;
  }

  return (double) s[i].et[0] + (double) s[i].et[1] * 1e-6;
}

double Timer_read_busy(int *clock)
{
  int i = *clock - 1;

  if (i < 0 || i >= NCLOCKS) {
    return 0.0;
  }

  return (double) s[i].ct[0] + (double) s[i].ct[1] * 1.0e-6;
}

double Timer_read_idle(int *clock)
{
  int i = *clock - 1;

  if (i < 0 || i >= NCLOCKS) {
    return 0.0;
  }

  return (double) s[i].it[0] + (double) s[i].it[1] * 1.0e-6;
}
