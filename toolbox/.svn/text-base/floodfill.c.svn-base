/**
 * floodfill.c: this file is part of the ALNSB project.
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

#include <toolbox/floodfill.h>

struct s_pair {
  int i;
  int j;
};
typedef struct s_pair s_pair_t;

static
void myfloodfill(int startRow, int startCol, int rows, int cols,
		 ALNSB_IMAGE_TYPE_BIN im[rows][cols],
		 ALNSB_IMAGE_TYPE_BIN newColor,
		 ALNSB_IMAGE_TYPE_BIN oldColor)
{
  if (im[startRow][startCol] != oldColor)
    return;
  s_pair_t* queue = (s_pair_t*) malloc (sizeof(s_pair_t) * rows*cols);
  queue[0].i = startRow;
  queue[0].j = startCol;
  unsigned int q_get_pos = 0;
  unsigned int q_push_pos = 1;
  int qs = rows * cols;
  while (q_get_pos != q_push_pos)
    {
      s_pair_t pRC;
      pRC.i = queue[q_get_pos % qs].i;
      pRC.j = queue[q_get_pos % qs].j;
      ++q_get_pos;
      if (im[pRC.i][pRC.j] == oldColor)
	{
	  int leftR = pRC.j;
	  int rightR = pRC.j;
	  while (leftR-1 >= 0 && im[pRC.i][leftR-1] == oldColor)
	    leftR--;
	  while (rightR+1 < cols && im[pRC.i][rightR] == oldColor)
	    rightR++;
	  int j;
	  for (j = leftR; j <= rightR; ++j)
	    {
	      im[pRC.i][j] = newColor;
	      if (pRC.i-1 >= 0 && im[pRC.i-1][j] == oldColor)
		{
		  queue[q_push_pos % qs].i = pRC.i-1;
		  queue[q_push_pos % qs].j = j;
		  ++q_push_pos;
		}
	      if (pRC.i+1 < cols && im[pRC.i+1][j] == oldColor)
		{
		  queue[q_push_pos % qs].i = pRC.i+1;
		  queue[q_push_pos % qs].j = j;
		  ++q_push_pos;
		}
	    }
	}
    }
  free (queue);
}


void alnsb_inplace_imfill_bin2d(int rows, int cols,
				ALNSB_IMAGE_TYPE_BIN im[rows][cols])
{
  ALNSB_IMAGE_TYPE_BIN newColor = 1;
  ALNSB_IMAGE_TYPE_BIN oldColor = 0;
  ALNSB_IMAGE_TYPE_BIN* imcopy =
    (ALNSB_IMAGE_TYPE_BIN*) malloc (sizeof(ALNSB_IMAGE_TYPE_BIN) * rows * cols);
  ALNSB_IMAGE_TYPE_BIN (*imcopy2d)[cols] =
    (ALNSB_IMAGE_TYPE_BIN (*)[cols])imcopy;
  unsigned int i, j;
  for (i = 0; i < rows; ++i)
    for (j = 0; j < cols; ++j)
      imcopy2d[i][j] = im[i][j];
  myfloodfill (0, 0, rows, cols, imcopy2d, newColor, oldColor);
  myfloodfill (0, cols-1, rows, cols, imcopy2d, newColor, oldColor);
  for (i = 0; i < rows; ++i)
    for (j = 0; j < cols; ++j)
      if (imcopy2d[i][j] != 1)
      	im[i][j] = 1;
  free (imcopy);
}
