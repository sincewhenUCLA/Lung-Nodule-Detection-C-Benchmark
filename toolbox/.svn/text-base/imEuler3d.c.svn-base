/**
 * imEuler3d.c: this file is part of the ALNSB project.
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

#include <toolbox/imEuler3d.h>

int
alnsb_imEuler3d_bin3d(ALNSB_IMAGE_TYPE_BIN* __ALNSB_RESTRICT_PTR im3d,
		      int dim0, int dim1, int dim2)
{
   int i, j, k;
   int v = 0, e1 = 0, e2 = 0, e3 = 0, f1 = 0, f2 = 0, f3 = 0, s = 0;
   ALNSB_IMAGE_TYPE_BIN (*img)[dim1][dim2] =
     (ALNSB_IMAGE_TYPE_BIN (*)[dim1][dim2])im3d;
   int chi;

   /* Count the number of vertices */ 
#pragma omp parallel for private(j,k) reduction(+:v)
   for (i=0; i<dim0; i++) 
      for (j=0; j<dim1; j++) 
	 for (k=0; k<dim2; k++) 
	     v += img[i][j][k];

#pragma omp parallel for private(j,k) reduction(+:e1)
   for (i=0; i<dim0-1; i++)
      for (j=0; j<dim1; j++)
         for (k=0; k<dim2; k++)
              e1 += (img[i][j][k] != 0 && img[i+1][j][k] != 0);


#pragma omp parallel for private(j,k) reduction(+:e2)
   for (i=0; i<dim0; i++) 
      for (j=0; j<dim1-1; j++)
	 for (k=0; k<dim2; k++)
              e2 += (img[i][j][k] != 0 && img[i][j+1][k] != 0);

#pragma omp parallel for private(j,k) reduction(+:e3)
   for (i=0; i<dim0; i++)
      for (j=0; j<dim1; j++)
	 for (k=0; k<dim2-1; k++)
	      e3 += (img[i][j][k] != 0 && img[i][j][k+1] != 0);

#pragma omp parallel for private(j,k) reduction(+:f1)
    for (i=0; i<dim0; i++)
      for (j=0; j<dim1-1; j++)
	 for (k=0; k<dim2-1; k++) 
	     f1 += (img[i][j][k] != 0 && img[i][j+1][k] != 0 && img[i][j][k+1] != 0 && img[i][j+1][k+1] != 0);
         
#pragma omp parallel for private(j,k) reduction(+:f2)
    for (i=0; i<dim0-1; i++)
      for (j=0; j<dim1; j++)
         for (k=0; k<dim2-1; k++)
            f2 += (img[i][j][k] != 0 && img[i+1][j][k] != 0 && img[i][j][k+1] != 0 && img[i+1][j][k+1] != 0);

#pragma omp parallel for private(j,k) reduction(+:f3)
    for (i=0; i<dim0-1; i++)
      for (j=0; j<dim1-1; j++)
         for (k=0; k<dim2; k++)
            f3 += (img[i][j][k] != 0 && img[i+1][j][k] != 0 && img[i][j+1][k] != 0 && img[i+1][j+1][k] != 0);

#pragma omp parallel for private(j,k) reduction(+:s)
   for (k=0; k<dim2-1; k++)
      for (j=0; j<dim1-1; j++)
         for (i=0; i<dim0-1; i++)
	    s += (img[i][j][k] != 0 && img[i][j+1][k] != 0 && img[i+1][j][k] != 0 && img[i+1][j+1][k] != 0 && img[i][j][k+1] != 0 && img[i][j+1][k+1] != 0 && img[i+1][j][k+1] != 0 && img[i+1][j+1][k+1] != 0);

   chi = v - (e1 + e2 + e3) + (f1 + f2 + f3) - s;
   return chi; 
}
