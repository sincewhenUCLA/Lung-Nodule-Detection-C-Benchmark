/**
 * rotate.c: this file is part of the ALNSB project.
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
#include <toolbox/rotate.h>

/**
 * Rotate a FP image.
 *
 */
void alnsb_rotate_slices_real_img (image3DReal* __ALNSB_RESTRICT_PTR img)
{
  unsigned int i, j, k;

  image3DReal* tmp = (image3DReal*) image3D_duplicate (img->image3D);
  ALNSB_IMRealTo3D(img,im3D);
  ALNSB_IMRealTo3D(tmp,tmp3d);
  
#pragma omp parallel for private(j,k)
  for (i = 0; i < img->slices; ++i)
    for (j = 0; j < img->rows; ++j)
      for (k = 0; k < img->cols; ++k)
	{
	  tmp3d[i][k][img->rows - j - 1] = im3D[i][j][k];
	  /* ALNSB_IMAGE_TYPE_REAL tmp = im3D[i][j][k]; */
	  /* im3D[i][j][k] = im3D[i][k][j]; */
	  /* im3D[i][k][j] = tmp; */
	}
  image3D_copy (img->image3D, tmp->image3D);
  image3D_free (tmp->image3D);
}
