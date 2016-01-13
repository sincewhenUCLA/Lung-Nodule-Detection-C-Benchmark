/**
 * step.h: this file is part of the ALNSB project.
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
#ifndef ALNSB_STEP_H
# define ALNSB_STEP_H

# include <utilities/types.h>
# include <utilities/images.h>


struct alnsb_step {
  // Input image info.
  size_t		num_reads;
  image3D**		read;
  size_t		num_writes;
  image3D**		write;
};
typedef struct alnsb_step s_alnsb_step_t;

#define ALNSB_STEP_ADD_INPUT_IMG 1
#define ALNSB_STEP_ADD_OUTPUT_IMG 2

extern
void alnsb_step_add_io (s_alnsb_step_t** s, image3D* img, int where);

extern
void alnsb_step_push_input (s_alnsb_step_t** s, image3D* img);

extern
void alnsb_step_push_output (s_alnsb_step_t** s, image3D* img);

extern
void alnsb_step_free (s_alnsb_step_t* s);

extern
s_alnsb_step_t* alnsb_step_alloc ();




#endif //!ALNSB_STEP_H
