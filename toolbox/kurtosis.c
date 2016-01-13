/**
 * kurtosis.c: this file is part of the ALNSB project.
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
#include <float.h>
#include <math.h>

#include <toolbox/kurtosis.h>

float
alnsb_kurtosis_real1d(ALNSB_IMAGE_TYPE_REAL* __ALNSB_RESTRICT_PTR x, int dim0)
{
   int i;
   float retval = 0;
   float mean = 0, m1 = 0.0, m2 = 0.0;
   /* Compute the mean */
   for (i=0; i<dim0; i++) 
     mean += x[i];
   mean = mean / dim0;
   
   for (i=0; i<dim0; i++) 
     m1 += pow ((x[i] - mean), 4);
   m1 = m1 / dim0;
   
   for (i=0; i<dim0; i++) 
     m2 += pow ((x[i] - mean), 2);
   m2 = m2 / dim0;
   
   retval = m1 / (pow (m2, 2));

   return retval;
}
