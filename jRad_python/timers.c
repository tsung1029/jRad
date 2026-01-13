/***************************************************************************************************
  Default timers (use POSIX timers)
***************************************************************************************************/

/*
  Use standard POSIX timers based on gettimeofday. (This is actually how
  LAM and OpenMPI implement MPI_Wtime). Minimum resolution is 1 microssecond.
  
  The actual resolution is measured by finding the minimum difference >0 
  between succsessive gettimeofday calls. 
*/
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <sys/resource.h>

#include <sys/time.h>
#include <stdint.h>

int64_t timer_ticks_
(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  
  return ((int64_t)tv.tv_sec)*1000000 + (int64_t)tv.tv_usec;
}

double timer_interval_seconds_
(int64_t *start, int64_t *end)
{
  return (*end - *start) * 1.0e-6;
}

double timer_cpu_seconds_
( void )
{
    struct timeval tv;
    double wtime;
    gettimeofday(&tv, NULL);
    wtime = tv.tv_sec;
    wtime += (double)tv.tv_usec * 1.0e-6;
    return wtime;
}

double timer_resolution_
( void )
{
  struct timeval tv1, tv2;
   
  gettimeofday(&tv1, NULL);
  do {
       gettimeofday(&tv2, NULL);
  } while (tv1.tv_usec == tv2.tv_usec);
  
  return (double)(tv2.tv_usec - tv1.tv_usec) * 1.0e-6;
}
