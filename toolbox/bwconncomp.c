/**
 * bwconncomp.c: this file is part of the ALNSB project.
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
/*
 * @filename: bwconncomp.c
 * Author:
 * Louis-Noel Pouchet <pouchet@cs.ucla.edu>. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#undef min
#undef max
#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

#include <toolbox/bwconncomp.h>

static
void* protected_malloc(unsigned int sz)
{
  unsigned nitem = sz / sizeof(char);
  void* ret = calloc (nitem, sizeof(char));
  if (ret == NULL)
    {
      fprintf (stderr, "[ERROR][bwconncomp] Memory exhausted\n");
      exit (1);
    }
  return ret;
}

/**
 * Emulates bwconncomp(im, 18) for a 3D image of size dim1 x dim2 x dim3.
 *
 * returns (through pass-by-reference on addresses of NULL pointers,
 * if NULL is passed this output is not generated):
 *
 * - components_coordinates, an array int[num_components+1][...] with,
      in each row 'i', the position of each pixel belonging to
      component 'i'
 * - components_size, an array int[num_components+1] with
      in each row 'i', the number of pixels in the component.
 * - num_components, the number of connected components.
 * - labeled_image, the integer image with the component id for each
 * non-zero pixel in the input image.
 *
 *
 */


static
void alnsb_bwconncomp_bin_safe (ALNSB_IMAGE_TYPE_BIN* in_data,
			       int dim1, int dim2, int dim3,
			       int*** components_coordinates,
			       int** components_size,
			       int* num_components,
			      int** labeled_image,
			       int use_full_neighb)
{
  /* printf ("start bwconncomp\n"); */
  int debug = 0;
  unsigned int i, j, k;
  unsigned int sz = dim1 * dim2 * dim3;

  int* out_data =
    (int*) protected_malloc (sizeof(ALNSB_IMAGE_TYPE_BIN) * sz);
  int (*out3d)[dim2][dim3] =
    (int (*)[dim2][dim3])out_data;
  ALNSB_IMAGE_TYPE_BIN (*im3d)[dim2][dim3] =
    (ALNSB_IMAGE_TYPE_BIN (*)[dim2][dim3])in_data;

  int* parts = (int*) protected_malloc (sizeof(int) * sz);
  int* relabel = (int*) protected_malloc (sizeof(int) * sz);

  for (i = 0; i < sz; ++i)
    {
      out_data[i] = 0;
      parts[i] = 0;
    }
  unsigned int label_cnt = 1;
  int not_converged;
  unsigned int ii, jj, kk;
  int BS = 32;
  int curr_comp = 1;
  int curr_label_pos = 1;
  int iter = 0;
  do
    {
      ++iter;
      if (debug)
	{
	  printf ("bwconncomp: iter (%d x %d x %d)\n", dim1, dim2, dim3);
	  printf ("curr_label_pos=%d, curr_comp=%d\n", curr_label_pos, curr_comp);
	}
      not_converged = 0;
      for (ii = 0; ii < dim1/BS + 1; ++ii)
      	for (jj = 0; jj < dim2/BS + 1; ++jj)
      	  for (kk = 0; kk < dim3/BS + 1; ++kk)
      	    for (i = ii*BS; i < min((ii+1)*BS,dim1); ++i)
      	      for (j = jj*BS; j < min((jj+1)*BS,dim2); ++j)
      		for (k = kk*BS; k < min((kk+1)*BS,dim3); ++k)
		  {
		    if (im3d[i][j][k] == 0)
		      continue;
		    int l, m, n;
		    int Halo = 1;
		    int lmin = -1;
		    int need_up = 0;
		    int lb1 = max(((int)i-Halo),((int)0));
		    int ub1 = min(((int)i+Halo),((int)dim1-1));
		    int lb2 = max((int)j - Halo, (int)0);
		    int ub2 = min((int)j + Halo,(int)dim2-1);
		    int lb3 = max((int)k - Halo, (int)0);
		    int ub3 = min((int)k + Halo,(int)dim3-1);
		    for (l = lb1; l <= ub1; ++l)
		      for (m = lb2; m <= ub2; ++m)
		    	for (n = lb3; n <= ub3; ++n)
		    	  {
			    if (! use_full_neighb)

			      {
				if (dim1 != 1)
				  {
				    // check for 18-connected, that is
				    // the pixel (l,m,n) is not a
				    // corner of the cube neighborhood
				    // of (i,j,k):

				    if (l == i - Halo &&
					m == j - Halo &&
					n == k - Halo)
				      continue;

				    if (l == i + Halo &&
					m == j - Halo &&
					n == k - Halo)
				      continue;

				    if (l == i - Halo &&
					m == j + Halo &&
					n == k - Halo)
				      continue;

				    if (l == i + Halo &&
					m == j + Halo &&
					n == k - Halo)
				      continue;

				    if (l == i - Halo &&
					m == j - Halo &&
					n == k + Halo)
				      continue;

				    if (l == i + Halo &&
					m == j - Halo &&
					n == k + Halo)
				      continue;

				    if (l == i - Halo &&
					m == j + Halo &&
					n == k + Halo)
				      continue;

				    if (l == i + Halo &&
					m == j + Halo &&
					n == k + Halo)
				      continue;

				  }

			      }

		    	    // find smallest value in neighb.
		    	    if (out3d[l][m][n] != 0)
		    	      {
		    		int val = out3d[l][m][n];
		    		assert(val > 0);
		    		if (lmin > val)
		    		  need_up = 1;
		    		if (lmin == -1 || lmin > val)
				  lmin = val;
		    	      }
		    	  }
		    // no neighbor, use a new class.
		    if (lmin == -1)
		      out3d[i][j][k] = curr_comp++;
		    else
		      {
			out3d[i][j][k] = lmin;
		    	if (need_up)
		    	  {
		    	    // all ptrs in neighb must be reset
			    for (l = lb1; l <= ub1; ++l)
			      for (m = lb2; m <= ub2; ++m)
				for (n = lb3; n <= ub3; ++n)
		    	    	  {
				    if (! use_full_neighb)

				      if (dim1 != 1)
					{
					  // check for 18-connected, that is
					  // the pixel (l,m,n) is not a
					  // corner of the cube neighborhood
					  // of (i,j,k):

					  if (l == i - Halo &&
					      m == j - Halo &&
					      n == k - Halo)
					    continue;

					  if (l == i + Halo &&
					      m == j - Halo &&
					      n == k - Halo)
					    continue;

					  if (l == i - Halo &&
					      m == j + Halo &&
					      n == k - Halo)
					    continue;

					  if (l == i + Halo &&
					      m == j + Halo &&
					      n == k - Halo)
					    continue;

					  if (l == i - Halo &&
					      m == j - Halo &&
					      n == k + Halo)
					    continue;

					  if (l == i + Halo &&
					      m == j - Halo &&
					      n == k + Halo)
					    continue;

					  if (l == i - Halo &&
					      m == j + Halo &&
					      n == k + Halo)
					    continue;

					  if (l == i + Halo &&
					      m == j + Halo &&
					      n == k + Halo)
					    continue;

					}

		    	    	    // find smallest value in neighb.
		    	    	    if (out3d[l][m][n] != 0)
				      {
					out3d[l][m][n] = lmin;
		    	    		not_converged++;
		    	    	      }
		    	    	  }
		    	  }
		      }
		  }
      if (debug)
	printf ("bwconn: done iter (%d changed)\n", not_converged);
    }
  while (not_converged);

  int num_labels = 0;
  int num_bg = 0;
  for (i = 0; i < sz; ++i)
    if (out_data[i])
      {
	if ((parts[out_data[i]])++ == 0)
	  ++num_labels;
      }
    else
      ++num_bg;

  *num_components = num_labels;
  if (debug)
    printf ("found %d connected components\n", num_labels);
  if (labeled_image == NULL && components_coordinates == NULL &&
      components_size == NULL)
    {
      free (out_data);
      free (parts);
      free (relabel);
      return;
    }

  int* num_pix = NULL;
  if (components_size || components_coordinates)
    {
      // 3. Scan the image, update the components, count how many pixel in
      // each component.
      num_pix = (int*) protected_malloc (sizeof(int) * (num_labels + 1));
      for (i = 0; i <= num_labels; ++i)
  	num_pix[i] = 0;
      int pos = 0;
      for (i = 0; i < sz; ++i)
	if (parts[i] != 0)
	  {
	    relabel[i] = pos;
	    num_pix[pos++] = parts[i];
	  }
      //      num_pix[0] = num_bg;
      if (components_size)
  	*components_size = num_pix;
    }

  if (components_coordinates)
    {
      // 4. Format the connected component output.
      int** ret_coord =
  	(int**) protected_malloc (sizeof(int*) * (num_labels + 1));
      int* ret_coord_pos =
  	(int*) protected_malloc (sizeof(int) * (num_labels + 1));
      for (i = 0; i <= num_labels; ++i)
  	{
  	  ret_coord[i] = (int*) protected_malloc (sizeof(int) * (num_pix[i] + 1));
  	  ret_coord_pos[i] = 0;
  	}

      // 5. Re-scan the image, to produce the output pixel list per component.
      for (i = 0; i < sz; ++i)
  	{
	  /* if (out_data[i] == 0) */
	  /*   ret_coord[0][ret_coord_pos[0]++] = i; */
	  /* else */
	  if (out_data[i] != 0)
	    {
	      int true_label = relabel[out_data[i]];
	      ret_coord[true_label][ret_coord_pos[true_label]++] = i;
	    }
  	}
      *components_coordinates = ret_coord;
      free (ret_coord_pos);
    }

  if (labeled_image)
    {
      int* ret_img =
  	(int*) protected_malloc (sizeof(int) * sz);
      for (i = 0; i < sz; ++i)
	if (out_data[i] != 0)
	  ret_img[i] = relabel[out_data[i]];
	else
	  ret_img[i] = 0;
      *labeled_image = ret_img;
    }
  /* printf ("found %d connected component\n", *num_components); */
  /* for (i = 0; i < *num_components; ++i) */
  /*   printf ("sz[%d]=%d\n", i, (*components_size)[i]); */

  free (parts);
  free (relabel);
  free(out_data);
  /* printf ("done bwconncomp\n"); */
}






static
void alnsb_bwconncomp_bin_fast (ALNSB_IMAGE_TYPE_BIN* in_data,
			       int dim1, int dim2, int dim3,
			       int*** components_coordinates,
			       int** components_size,
			       int* num_components,
			       int** labeled_image,
			       int use_full_neighb)
{
  /// FIXME: port from dev branch.



}


void alnsb_bwconncomp_bin (ALNSB_IMAGE_TYPE_BIN* in_data,
			  int dim1, int dim2, int dim3,
			  int*** components_coordinates,
			  int** components_size,
			  int* num_components,
			  int** labeled_image,
			  int use_full_neighb)
{
  /* alnsb_bwconncomp_bin_fast (in_data, dim1, dim2, dim3, components_coordinates, */
  /* 			    components_size, num_components, labeled_image, */
  /* 			    use_full_neighb); */
  alnsb_bwconncomp_bin_safe (in_data, dim1, dim2, dim3, components_coordinates,
  			    components_size, num_components, labeled_image,
  			    use_full_neighb);
}
