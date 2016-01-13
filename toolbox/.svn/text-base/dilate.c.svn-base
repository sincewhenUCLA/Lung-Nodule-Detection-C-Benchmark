/**
 * dilate.c: this file is part of the ALNSB project.
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

#include <toolbox/dilate.h>
#include <toolbox/structuring_elements.h>


static
void dilate_im (ALNSB_IMAGE_TYPE_BIN* im, ALNSB_IMAGE_TYPE_BIN* out_im, int dim1, int dim2, int stw)
{
  ALNSB_IMAGE_TYPE_BIN (*im2d)[dim2] = (ALNSB_IMAGE_TYPE_BIN (*)[dim2])(im);
  ALNSB_IMAGE_TYPE_BIN (*out2d)[dim2] = (ALNSB_IMAGE_TYPE_BIN (*)[dim2])(out_im);
  int i, j;
  int halo = structuringElements[stw].halo;
  int width = 2 * halo + 1;
  int SE[width][width];
  for (i = 0; i < width; ++i) 
    for (j = 0; j < width; ++j)
      SE[i][j] = structuringElements[stw].se[i * (width) + j];

#pragma omp parallel for private(j) firstprivate(SE,halo,dim1,dim2,width)
  for (i = halo; i < dim1 - halo; ++i)
    for (j = halo; j < dim2 - halo; ++j)
      {
	int ii, jj;
	int sum = 0;
	for (ii = 0; ii < width; ++ii)
	  {
	    for (jj = 0; jj < width; ++jj)
	      {
		if (SE[ii][jj] == 0)
		  continue;
		if (im2d[i+ii-halo][j+jj-halo] != 0)
		  break;
	      }
	    if (jj < width)
	      break;
	  }
	if (ii < width)
	  out2d[i][j] = 1;
	else
	  out2d[i][j] = 0;
      }
}



void alnsb_inplace_dilate_bin2d(ALNSB_IMAGE_TYPE_BIN* im2d,
				int rows, int cols,
				int structelt)
{
  ALNSB_IMAGE_TYPE_BIN* out =
    (ALNSB_IMAGE_TYPE_BIN*) malloc (sizeof(ALNSB_IMAGE_TYPE_BIN) * rows * cols);
  dilate_im (im2d, out, rows, cols, structelt);
  int i;
  for (i = 0; i < rows * cols; ++i)
    im2d[i] = out[i];
  free (out);
}
