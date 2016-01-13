/**
 * preselection_step.c: this file is part of the ALNSB project.
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
#include <math.h>
#include <assert.h>
#include <stages/preselection/preselection_step.h>
#include <toolbox/bwconncomp.h>
#include <toolbox/imPerimeter.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif
#ifdef max
# undef max
#endif
#ifdef min
# undef min
#endif
#define min(a,b) (a < b ? a : b)
#define max(a,b) (a > b ? a : b)



void preselection_cpu (s_alnsb_environment_t* __ALNSB_RESTRICT_PTR env,
		       image3DBin* __ALNSB_RESTRICT_PTR input,
		       image3DBin** __ALNSB_RESTRICT_PTR output)
{
  // Allocate output, and copy input image to it.
  *output = image3DBin_alloc (input->slices, input->rows, input->cols);
  ALNSB_IMBinTo1D(*output,out_img);
  // View input as 1D image.
  ALNSB_IMBinTo1D(input,im_in);

  // Algorithm parameters.
  float diameTMin = env->preselection_diameterMin;
  float diameTMax = env->preselection_diameterMax;
  float elongationTMax = env->preselection_elongationMax;
  float circulTMin = env->preselection_circulMin;

  // Local variables.
  int heights = input->slices;
  int rows = input->rows;
  int cols = input->cols;
  float space_x = env->scanner_pixel_spacing_x_mm;
  float space_y = env->scanner_pixel_spacing_x_mm;
  float space_z = env->scanner_slice_thickness_mm;
  int sz = heights * rows * cols;
  int i;
  float areaTMin = pow((diameTMin/2.0),2) * M_PI;
  float areaTMax = pow((diameTMax/2.0),2) * M_PI;
  float volumTMin = 3*pow((diameTMin/2.0),3) * M_PI/4.0;
  float volumTMax = 3*pow((diameTMax/2.0),3) * M_PI/4.0;

  int nodules_count = 0;
  int** comp_coordinates = NULL;
  int* comp_sz = NULL;
  int nbcomp = 0;
  alnsb_bwconncomp_bin (im_in, heights, rows, cols,
			&comp_coordinates, &comp_sz, &nbcomp, NULL, 0);

  if (env->verbose_level > 0)
    printf ("[INFO] Preselection on %d candidate objects\n", nbcomp);
#pragma omp parallel for
  for (i = 0; i < nbcomp; ++i)
    {
      if (comp_sz[i] < 5)
	continue;
      int j, k;
      // get min/max coord in each dimension.
      int min_x = sz, min_y = sz, min_z = sz;
      int max_x = -1, max_y = -1, max_z = -1;
      for (j = 0; j < comp_sz[i]; ++j)
	{
	  // convert coord. to 3D map.
	  int z = comp_coordinates[i][j] / (rows * cols);
	  int x = (comp_coordinates[i][j] % (rows * cols)) / cols;
	  int y = (comp_coordinates[i][j] % (rows * cols)) % cols;
	  min_z = z < min_z ? z : min_z;
	  min_x = x < min_x ? x : min_x;
	  min_y = y < min_y ? y : min_y;
	  max_z = z > max_z ? z : max_z;
	  max_x = x > max_x ? x : max_x;
	  max_y = y > max_y ? y : max_y;
	}
      float xLength = (max_x - min_x + 1) * space_x;
      float yLength = (max_y - min_y + 1) * space_y;
      float zLength = (max_z - min_z + 1) * space_z;
      float diameter = max(xLength, yLength);
      diameter = max(diameter, zLength);
      float minSz = min(xLength, yLength);
      minSz = min(minSz,zLength);
      float elongation = diameter/minSz;
      float volume = comp_sz[i] * space_x * space_y * space_z;
      // LNP: Compute the area of the 2D slice centered along the z
      // axis for the component.
      int z_idx = (max_z + min_z) / 2;

      if (diameter < diameTMin)
	continue;
      if (diameter > diameTMax)
	continue;
      if (elongation > elongationTMax)
	continue;
      if (volume > volumTMax)
	continue;
      if (volume < volumTMin)
	continue;

      ALNSB_IMAGE_TYPE_BIN* im_slice =
	(ALNSB_IMAGE_TYPE_BIN*) malloc (sizeof(ALNSB_IMAGE_TYPE_BIN) *rows*cols);
      for (j = 0; j < rows * cols; ++j)
	im_slice[j] = 0;
      int has_px = 0;
      for (j = 0; j < comp_sz[i]; ++j)
	{
	  int coord = comp_coordinates[i][j];
	  if (coord / (rows * cols) == z_idx)
	    {
	      int planepos = coord % (rows * cols);
	      im_slice[planepos] = 1;
	      has_px = 1;
	    }
	}
      if (! has_px)
	{
	  assert(0);
	  has_px = 0;
	  // safety net, use min_z for the base slice.
	  for (j = 0; j < comp_sz[i]; ++j)
	    {
	      int coord = comp_coordinates[i][j];
	      if (coord / (rows * cols) == min_z)
		{
		  int planepos = coord % (rows * cols);
		  im_slice[planepos] = 1;
		  has_px = 1;
		}
	    }
	}
      assert(has_px);
      int** comp_coordinatesX = NULL;
      int* comp_szX = NULL;
      int nbcompX = 0;
      alnsb_bwconncomp_bin (im_slice, 1, rows, cols,
			    &comp_coordinatesX, &comp_szX, &nbcompX, NULL, 1);
      assert (nbcompX > 0);
      int max_sz = -1; int max_id;
      for (j = 0; j < nbcompX; ++j)
	if (max_sz < comp_szX[j])
	  {
	    max_sz = comp_szX[j];
	    max_id = j;
	  }
      assert(max_id >= 0);
      for (j = 0; j < nbcompX; ++j)
	{
	  if (j == max_id)
	    continue;
	  for (k = 0; k < comp_szX[j]; ++k)
	    im_slice[comp_coordinatesX[j][k]] = 0;
	  free (comp_coordinatesX[j]);
	}
      float area = (float) comp_szX[max_id];
      free (comp_coordinatesX[max_id]);
      free (comp_szX);
      free (comp_coordinatesX);
      float perimeter = alnsb_imPerimeter_bin2d (im_slice, rows, cols);
      float roundDegree = (4 * M_PI * area) / pow (perimeter, 2);
      free (im_slice);
      if (area > areaTMax || area < areaTMin)
	continue;
      if (roundDegree < circulTMin)
      	continue;

      // Good candidate nodule.
      for (j = 0; j < comp_sz[i]; ++j)
      	out_img[comp_coordinates[i][j]] = 1;
      #pragma atomic
      ++nodules_count;
    }

  if (env->verbose_level > 0)
    printf ("[INFO] Preselection retained %d candidate nodules\n",
	    nodules_count);

  for (i = 0; i < nbcomp; ++i)
    free (comp_coordinates[i]);
  free (comp_coordinates);
  free (comp_sz);
}
