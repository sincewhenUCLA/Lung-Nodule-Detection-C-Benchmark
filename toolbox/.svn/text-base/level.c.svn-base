/**
 * level.c: this file is part of the ALNSB project.
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
#include <toolbox/level.h>


/**
 * Level a FP image.
 *
 * img(img > upperBand) = upperBand
 * img(img < lowerBand) = lowerBand
 *
 */
void alnsb_level_real_img (image3DReal* __ALNSB_RESTRICT_PTR img,
			   ALNSB_IMAGE_TYPE_REAL lowerBand,
			   ALNSB_IMAGE_TYPE_REAL upperBand)
{
  unsigned int i;

  ALNSB_IMRealTo1D(img, im1D);

#pragma omp parallel for
  for (i = 0; i < img->num_pixels; ++i)
    {
      im1D[i] = im1D[i] > upperBand ? upperBand : im1D[i];
      im1D[i] = im1D[i] < lowerBand ? lowerBand : im1D[i];
    }
}
