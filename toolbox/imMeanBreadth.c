/**
 * imMeanBreadth.c: this file is part of the ALNSB project.
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

#include <toolbox/imMeanBreadth.h>

float
alnsb_imMeanBreadth_bin3d(ALNSB_IMAGE_TYPE_BIN* __ALNSB_RESTRICT_PTR im3d,
			  float *xyzSpace, int dim0, int dim1, int dim2)
{
   int i, j, k;
   int nv = 0, ne1 = 0, ne2 = 0, ne3 = 0, nf1 = 0, nf2 = 0, nf3 = 0;
   int b1, b2, b3;
   float breadth, a1, a2, a3;

   ALNSB_IMAGE_TYPE_BIN (*img)[dim1][dim2] =
     (ALNSB_IMAGE_TYPE_BIN (*)[dim1][dim2])im3d;
   float d1 = xyzSpace[0];
   float d2 = xyzSpace[1];
   float d3 = xyzSpace[2];
 
#pragma omp parallel for reduction(+:nv)
   for (i=0; i<dim0; i++) 
      for (j=0; j<dim1; j++) 
	 for (k=0; k<dim2; k++) 
	     nv += img[i][j][k] != 0;

#pragma omp parallel for reduction(+:ne2)
   for (i=0; i<dim0; i++) 
      for (j=0; j<dim1-1; j++)
	 for (k=0; k<dim2; k++)
              ne2 += (img[i][j][k] != 0 && img[i][j+1][k] != 0);

#pragma omp parallel for reduction(+:ne1)
   for (i=0; i<dim0-1; i++)
      for (j=0; j<dim1; j++)
	 for (k=0; k<dim2; k++)
	      ne1 += (img[i][j][k] != 0 && img[i+1][j][k] != 0);

#pragma omp parallel for reduction(+:ne3)
   for (i=0; i<dim0; i++)
      for (j=0; j<dim1; j++)
	 for (k=0; k<dim2-1; k++)
	      ne3 += (img[i][j][k] != 0 && img[i][j][k+1] != 0);

#pragma omp parallel for reduction(+:nf1)
    for (i=0; i<dim0; i++)
      for (j=0; j<dim1-1; j++)
	 for (k=0; k<dim2-1; k++) 
	     nf1 += (img[i][j][k] != 0 && img[i][j+1][k] != 0 && img[i][j][k+1] != 0 && img[i][j+1][k+1] != 0);
         
#pragma omp parallel for reduction(+:nf2)
    for (i=0; i<dim0-1; i++)
      for (j=0; j<dim1; j++)
         for (k=0; k<dim2-1; k++)
            nf2 += (img[i][j][k] != 0 && img[i+1][j][k] != 0 && img[i][j][k+1] != 0 && img[i+1][j][k+1] != 0);

#pragma omp parallel for reduction(+:nf3)
    for (i=0; i<dim0-1; i++)
      for (j=0; j<dim1-1; j++)
         for (k=0; k<dim2; k++)
            nf3 += (img[i][j][k] != 0 && img[i+1][j][k] != 0 && img[i][j+1][k] != 0 && img[i+1][j+1][k] != 0);

   b1 = nv - (ne2 + ne3) + nf1;
   b2 = nv - (ne1 + ne3) + nf2;
   b3 = nv - (ne1 + ne2) + nf3;

   a1 = (d1 * d2 * d3) / (d2 * d3);
   a2 = (d1 * d2 * d3) / (d1 * d3);
   a3 = (d1 * d2 * d3) / (d1 * d2);
  
   breadth = (b1 * a1 + b2 * a2 + b3 * a3)/3;
   return breadth; 
}
