/**
 * scale.c: this file is part of the ALNSB project.
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
#include <toolbox/scale.h>


/**
 * Scale a FP image.
 *
 * img = (img-min(img(:)))./(max(img(:))-min(img(:)))
 *
 */
void alnsb_scale_real_img (image3DReal* __ALNSB_RESTRICT_PTR img)
{
  unsigned int i;

  ALNSB_IMRealTo1D(img, im1D);

  ALNSB_IMAGE_TYPE_REAL min_im = im1D[0];
  ALNSB_IMAGE_TYPE_REAL max_im = im1D[0];
  for (i = 1; i < img->num_pixels; ++i)
    {
      ALNSB_IMAGE_TYPE_REAL px = im1D[i];
      min_im = px < min_im ? px : min_im;
      max_im = px > max_im ? px : max_im;
    }

#pragma omp parallel for
  for (i = 0; i < img->num_pixels; ++i)
    im1D[i] = (im1D[i] - min_im) / (max_im - min_im);
}
