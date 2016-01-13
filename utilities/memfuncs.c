/**
 * memfuncs: this file is part of the ALNSB project.
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

#include <utilities/common.h>


/**
 *
 *
 */
void* alnsb_calloc (size_t data_sz, size_t num_elt)
{
  void* data = calloc (num_elt, data_sz);
  if (data == NULL)
    {
      fprintf (stderr, "[%s][ERROR] Memory exhausted\n", ALNSB_PROJECT_STR);
      exit (1);
    }

  return data;
}
