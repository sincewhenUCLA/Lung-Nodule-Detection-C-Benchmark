/**
 * stdev.c: this file is part of the ALNSB project.
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
#include <math.h>
#include <toolbox/stdev.h>

float
alnsb_stdev_real1d(ALNSB_IMAGE_TYPE_REAL* __ALNSB_RESTRICT_PTR x, int dim0)
{
  int i;
  float retval = 0;
  float mean = 0.0, m1 = 0.0;
  /* Compute the mean */
  for (i=0; i<dim0; i++) 
    mean += x[i];
  mean = mean / (float)dim0;

  for (i=0; i<dim0; i++) 
    m1 += ((x[i] - mean)*(x[i] - mean));
  
  retval = sqrt (m1/(float)(dim0-1));
  retval = round (retval*100000)/100000; 
  return retval;
}
