/**
 * timer.c: this file is part of the ALNSB project.
 *
 * ALNSB: the Adaptive Lung Nodule Screening Benchmark
 *
 * Copyright (C) 2014,2015 University of California Los Angeles
 *
 * This program can be redistributed and/or modified under the terms
 * of the license specified in the LICENSE.txt file at the root of the
 * project.
 *
 * Contact: Alex Bui <buia@mii.ucla.edu>
 *
 */
/**
 * Written by: Shiwen Shen, Prashant Rawat, Louis-Noel Pouchet and William Hsu
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#ifdef _OPENMP
# include <omp.h>
#endif

#include <utilities/timer.h>
#include <utilities/common.h>


/////////////////////////////////////////////////////////////
////////// Timer, extracted from PolyBench/C 3.2. ///////////
/////////////////////////////////////////////////////////////
static double polybench_t_start;
static double polybench_t_end;

static double rtclock ()
{
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, NULL);
    if (stat != 0)
      printf ("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

static
void polybench_timer_start ()
{
  polybench_t_start = rtclock ();
}

static
void polybench_timer_stop ()
{
  polybench_t_end = rtclock ();
}

static
void polybench_timer_print (FILE* f)
{
  fprintf (f, "%0.3f seconds", polybench_t_end - polybench_t_start);
}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


void alnsb_timer_start ()
{
  polybench_timer_start ();
}


void alnsb_timer_stop ()
{
  polybench_timer_stop ();
}


void alnsb_timer_print (FILE* f, char* message)
{
  fprintf (f, "[%s][%s] ", ALNSB_PROJECT_STR, message);
  polybench_timer_print (f);
  fprintf (f, "\n");
}
