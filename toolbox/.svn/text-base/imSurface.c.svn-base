/**
 * imSurface.c: this file is part of the ALNSB project.
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

#include <toolbox/imSurface.h>

float
alnsb_imSurface_bin3d(ALNSB_IMAGE_TYPE_BIN* __ALNSB_RESTRICT_PTR im3d,
		      float *xyzSpace, int dim0, int dim1, int dim2)
{
   int i, j, k;
   int nv = 0, n1, n2, n3, n1t = 0, n2t = 0, n3t = 0;
   float surf, n1byd1, n2byd2, n3byd3;
   ALNSB_IMAGE_TYPE_BIN (*img)[dim1][dim2] =
     (ALNSB_IMAGE_TYPE_BIN (*)[dim1][dim2])im3d;
   float d1 = xyzSpace[0];
   float d2 = xyzSpace[1];
   float d3 = xyzSpace[2];

#pragma omp parallel for private(j,k) reduction(+:nv)
   for (i=0; i<dim0; i++)
      for (j=0; j<dim1; j++)
	 for (k=0; k<dim2; k++)
	     nv += img[i][j][k] != 0;

#pragma omp parallel for private(j,k) reduction(+:n1t)
   for (i=0; i<dim0; i++)
      for (j=0; j<dim1-1; j++)
	 for (k=0; k<dim2; k++)
              n1t += (img[i][j][k] != 0 && img[i][j+1][k] != 0);

#pragma omp parallel for private(j,k) reduction(+:n2t)
   for (i=0; i<dim0-1; i++)
      for (j=0; j<dim1; j++)
	 for (k=0; k<dim2; k++)
	      n2t += (img[i][j][k] != 0 && img[i+1][j][k] != 0);

#pragma omp parallel for private(j,k) reduction(+:n3t)
   for (i=0; i<dim0; i++)
      for (j=0; j<dim1; j++)
	 for (k=0; k<dim2-1; k++)
	      n3t += (img[i][j][k] != 0 & img[i][j][k+1] != 0);

   n1 = nv - n1t;
   n2 = nv - n2t;
   n3 = nv - n3t;

   n1byd1 = (float)n1/d1;
   n2byd2 = (float)n2/d2;
   n3byd3 = (float)n3/d3;

   surf = (4 * (n1byd1 + n2byd2 + n3byd3) * (d1 * d2 * d3))/3;
   return surf;
}
