/**
 * preparation_step.c: this file is part of the ALNSB project.
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
#include <stages/preparation/preparation_step.h>
#include <toolbox/rotate.h>
#include <toolbox/level.h>
#include <toolbox/scale.h>


void preparation_cpu (s_alnsb_environment_t* __ALNSB_RESTRICT_PTR env,
		      image3DReal* __ALNSB_RESTRICT_PTR input,
		      image3DReal** __ALNSB_RESTRICT_PTR output)
{
  if (env->verbose_level > 1)
    printf ("[DEBUG] Start pass preparation-cpu\n");

  // Allocate output, and copy input image to it.
  *output = (image3DReal*) image3D_duplicate (input->image3D);
  // Perform rotation 90 degrees + leveling + scaling.
  if (env->transpose_input)
    alnsb_rotate_slices_real_img (*output);
  alnsb_level_real_img (*output, env->thresold_lowerBand,
			env->thresold_upperBand);
  alnsb_scale_real_img (*output);

  if (env->verbose_level > 1)
    printf ("[DEBUG] Done pass preparation-cpu\n");
}
