/**
 * segmentationMask_step.c: this file is part of the ALNSB project.
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
#include <assert.h>


#include <stages/segmentationMask/segmentationMask_step.h>
#include <toolbox/floodfill.h>
#include <toolbox/dilate.h>
#include <toolbox/erode.h>
#include <toolbox/structuring_elements.h>
#include <toolbox/bwconncomp.h>

static
void get_two_largest_comp (int* comp_sz, int nb_comp, int* idx1, int* idx2)
{
  int m1 = 0;
  int id1 = -1;
  int id2 = -1;
  int i;
  // Get the largest.
  for (i = 0; i < nb_comp; ++i)
    {
      if (m1 < comp_sz[i])
	{
	  m1 = comp_sz[i];
	  id1 = i;
	}
    }
  // Get the 2nd largest.
  m1 = 0;
  for (i = 0; i < nb_comp; ++i)
    {
      if (i == id1)
	continue;
      if (m1 < comp_sz[i])
	{
	  m1 = comp_sz[i];
	  id2 = i;
	}
    }
  *idx1 = id1;
  *idx2 = id2;
}


void segmentationMask_cpu (s_alnsb_environment_t* __ALNSB_RESTRICT_PTR env,
			   image3DReal* __ALNSB_RESTRICT_PTR input,
			   image3DBin** __ALNSB_RESTRICT_PTR output)
{
  // Allocate output data.
  *output = image3DBin_alloc (input->slices, input->rows, input->cols);
  // 1D view of output data.
  ALNSB_IMBinTo1D(*output, output_image1d);

  // 1D view of input data.
  ALNSB_IMRealTo1D(input, seg_image1d);

  // Local variables.
  unsigned int rows = input->rows;
  unsigned int cols = input->cols;
  unsigned int heights = input->slices;
  unsigned int sz = rows * cols * heights;
  unsigned int i, j, k;

  // Local images.
  // Allocate final mask, use 1D view.
  ALNSB_BIN1D_alloc(output_mask1d, input->slices, input->rows, input->cols);

  // output_mask1d = binary(input)
  for (i = 0; i < sz; ++i)
    output_mask1d[i] = (seg_image1d[i] != 0);

  // floodfill2d(output_mask) on all 2D slices
  ALNSB_IMAGE_TYPE_BIN* ptr1 = output_mask1d;
  for (i = 0; i < heights; ++i)
    {
      ALNSB_2Dslice_from_Bin1D(slice2d, ptr1, rows, cols);
      alnsb_inplace_imfill_bin2d (rows, cols, slice2d);
      ptr1 += rows * cols;
    }

  // output_mask1d = seg_image1d && output_mask1d
  for (i = 0; i < sz; ++i)
    output_mask1d[i] = output_mask1d[i] == 0 ? 0 : seg_image1d[i] == 0;

  // Collect all 3D objects.
  int** comp_coordinates = NULL;
  int* comp_sz = NULL;
  int nbcomp = 0;
  int* mask = NULL;
  alnsb_bwconncomp_bin (output_mask1d, heights, rows, cols,
		       &comp_coordinates, &comp_sz, &nbcomp, &mask, 1);
  
  // Keep only the two largest 3D objects in the mask.
  int idx1, idx2;
  get_two_largest_comp (comp_sz, nbcomp, &idx1, &idx2);
  for (i = 0; i < sz; ++i)
    output_mask1d[i] = 0;
  if (idx1 >= 0)
    for (i = 0; i < comp_sz[idx1]; ++i)
      output_mask1d[comp_coordinates[idx1][i]] = 1;
  if (idx2 >= 0)
    for (i = 0; i < comp_sz[idx2]; ++i)
      output_mask1d[comp_coordinates[idx2][i]] = 1;
  for (i = 0; i < nbcomp; ++i)
    free (comp_coordinates[i]);
  free (comp_coordinates);
  free (comp_sz);


  // Proceed each 2D slice of the mask independently, and perform:
  //   - floodfill2d
  //   - closing
  //   - erode
  //   - take 2 largest objects, keep 1st, and 2nd if sz(2nd) > .3 * sz(1st)

  /* #pragma omp parallel for */
  for (i = 0; i < heights; ++i)
    {
      ALNSB_IMAGE_TYPE_BIN* ptrval = output_mask1d;
      ptrval += (i * rows * cols);
      int** comp_coordinatesL = NULL;
      int* comp_szL = NULL;
      int nbcompL = 0;
      int* maskL = NULL;
      ALNSB_2Dslice_from_Bin1D(im2d, ptrval, rows, cols);
      alnsb_inplace_imfill_bin2d (rows, cols, im2d);
      alnsb_inplace_dilate_bin2d (ptrval, rows, cols, SE_2D_diamond_5);
      alnsb_inplace_erode_bin2d (ptrval, rows, cols, SE_2D_diamond_5);
      alnsb_inplace_erode_bin2d (ptrval, rows, cols, SE_2D_diamond_2);

      alnsb_bwconncomp_bin (ptrval, 1, rows, cols,
			    &comp_coordinatesL, &comp_szL, &nbcompL, NULL, 1);
      if (nbcomp > 0)
      	{
      	  int idx1, idx2;
      	  get_two_largest_comp (comp_szL, nbcompL, &idx1, &idx2);
      	  float total1 = (float)comp_szL[idx1] / ((float)rows*cols);
      	  for (j = 0; j < rows * cols; ++j)
      	    output_mask1d[i*rows*cols + j] = 0;
      	  if (idx1 >= 0)
      	    for (j = 0; j < comp_szL[idx1]; ++j)
      	      output_mask1d[i*rows*cols + comp_coordinatesL[idx1][j]] = 1;
      	  if (idx2 >= 0)
      	    {
      	      float total2 = (float)comp_szL[idx2] / ((float)rows*cols);
      	      if ((total1 / total2) < 3)
      		for (j = 0; j < comp_szL[idx2]; ++j)
      		  output_mask1d[i*rows*cols + comp_coordinatesL[idx2][j]] = 1;
      	    }
      	}
      for (j = 0; j < nbcompL; ++j)
      	free (comp_coordinatesL[j]);
      free (comp_coordinatesL);
      free (comp_szL);
    }

  // Apply mask to input to form output binary image.
  for (i = 0; i < sz; ++i)
    if (output_mask1d[i] != 0 && seg_image1d[i] != 0)
      output_image1d[i] = 1;
    else
      output_image1d[i] = 0;
  // Free temporary mask.
  ALNSB_BIN1D_free(output_mask1d);

  // Remove first and last 4 slices, if asked.
  if (env->segmentationMask_skip_4_slices)
    {
      unsigned int slice4 = rows * cols * 4;
      for (i = 0; i < slice4; ++i)
  	output_image1d[i] = 0;
      for (i = sz - slice4; i < sz; ++i)
  	output_image1d[i] = 0;
    }
}
