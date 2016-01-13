/**
 * step.c: this file is part of the ALNSB project.
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
#include <stdlib.h>
#include <stdio.h>

#include <utilities/step.h>
#include <utilities/memfuncs.h>
#include <utilities/images.h>

s_alnsb_step_t* alnsb_step_alloc ()
{
  s_alnsb_step_t* ret = (s_alnsb_step_t*) malloc (sizeof(s_alnsb_step_t));
  if (! ret)
    {
      fprintf (stderr, "[ERROR] Memory exhausted\n");
      exit (1);
    }
  ret->num_reads = ret->num_writes = 0;
  ret->read = ret->write = NULL;

  return ret;
}

void alnsb_step_add_io (s_alnsb_step_t** s, image3D* img, int where)
{
  if (s == NULL)
    exit (1);
  if (*s == NULL)
    *s = alnsb_step_alloc ();
  image3D*** ptr;
  size_t* sz;
  if (where == ALNSB_STEP_ADD_INPUT_IMG)
    {
      ptr = &((*s)->read);
      sz = &((*s)->num_reads);
      if (img) img->ref_counter++;
    }
  else if (where == ALNSB_STEP_ADD_OUTPUT_IMG)
    {
      ptr = &((*s)->write);
      sz = &((*s)->num_writes);
    }
  (*(sz))++;
  *ptr = (image3D**) realloc (*ptr, *sz * sizeof(image3D*));
  (*ptr)[*sz - 1] = img;
}

void alnsb_step_push_input (s_alnsb_step_t** s, image3D* img)
{
  alnsb_step_add_io (s, img, ALNSB_STEP_ADD_INPUT_IMG);
}


void alnsb_step_push_output (s_alnsb_step_t** s, image3D* img)
{
  alnsb_step_add_io (s, img, ALNSB_STEP_ADD_OUTPUT_IMG);
}


void alnsb_step_free (s_alnsb_step_t* s)
{
  if (! s)
    return;
  int i;
  for (i = 0; i < s->num_reads; ++i)
    if (s->read[i])
      image3D_free (s->read[i]);
  free (s->read);
  for (i = 0; i < s->num_writes; ++i)
    if (s->write[i])
      image3D_free (s->write[i]);
  free (s->write);
  free (s);
}
