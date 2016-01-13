/**
 * imPerimeter.c: this file is part of the ALNSB project.
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
#include <toolbox/imPerimeter.h>

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

float
alnsb_imPerimeter_bin2d(ALNSB_IMAGE_TYPE_BIN* __ALNSB_RESTRICT_PTR im2d,
			int rows, int cols)
{
  /* ALNSB_IMAGE_TYPE_BIN (*img)[cols] = (ALNSB_IMAGE_TYPE_BIN (*)[cols])im2d; */

  /* unsigned int perim = 0; */
  /* unsigned int i, j; */
  /* for (i=1; i<rows-1; i++) */
  /*   for (j=1; j<cols-1; j++) */
  /*     if (img[i][j] != 0 && */
  /* 	  (img[i-1][j] == 0 || */
  /* 	   img[i][j-1] == 0 || */
  /* 	   img[i+1][j] == 0 || */
  /* 	   img[i][j+1] == 0)) */
  /* 	++perim; */
  
  /*  return perim; */

  /// LNP: imPerimeter code, edited to comply with
  /// http://www.mathworks.com/matlabcentral/fileexchange/33690-geometric-measures-in-2d-3d-images/content/imMinkowski/imPerimeterEstimate.m
  /// 
  /// It uses approximation while computing the exact perimeter (as
  /// above) is possible, so I don't know why we use this...
   
   int i, j, k;
   int nv = 0, n1, n2, n3, n4, n1t = 0, n2t = 0, n3t = 0, n4t = 0;
  ALNSB_IMAGE_TYPE_BIN (*img)[cols] = (ALNSB_IMAGE_TYPE_BIN (*)[cols])im2d;
   double perim, d1 = 1, d2 = 1;
   double d12 = sqrt (pow(d1,2) + pow(d2,2));
   double vol = d1 * d2;

#pragma omp parallel for private(i, j) reduction(+: nv, n1t, n2t, n3t, n4t)
   for (i=1; i<rows; i++)
	   for (j=1; j<cols; j++)  {
	     n1t += ((img[i-1][j] != 0) != (img[i][j] != 0));
	     n2t += ((img[i][j-1] != 0) != (img[i][j] != 0));
	     n3t += ((img[i-1][j-1] != 0) != (img[i][j] != 0));
	     n4t += ((img[i-1][j] != 0) != (img[i][j-1] != 0));
	   }

   perim = ((n1t*d1 + n2t*d2 + ((double)(n3t+n4t))/d12) * M_PI) / 8;

   return perim;
}
