/**
 * main.c: this file is part of the ALNSB project.
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

#include <utilities/environment.h>
#include <driver/pipeline.h>

extern int alnsb_getopts (int argc, char** argv, s_alnsb_environment_t* env);


int main (int argc, char** argv)
{
  s_alnsb_environment_t* env = alnsb_environment_malloc ();
  alnsb_getopts (argc, argv, env);
  
  alnsb_pipeline (env);

  alnsb_environment_free (env);

  return 0;
}
