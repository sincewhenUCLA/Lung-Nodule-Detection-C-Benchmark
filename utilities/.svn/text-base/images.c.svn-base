/**
 * images.c: this file is part of the ALNSB project.
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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <utilities/images.h>
#include <utilities/memfuncs.h>
#include <utilities/types.h>

/**
 *
 * Data layout: a 3D image is allocated as a 1-dimensional vector of
 * data, using the macro ALNSB_ALLOC3D(type, size1, size2, size3)
 * which returns a type* pointer.
 *
 * This pointer can be viewed as a 3D image im[i][j][k] using the
 * macro ALNSB_PTRTOARRAY3D(ptrname, 3darrname, type, size1, size2, size3)
 * to declare a new variable 3darrname which is a C99 3D array (e.g.,
 * im[i][j][k]) from a pointer ptrname (e.g., as returned by ALNSBD_ALLOC3D)
 * such that:
 *   - along the first dimension (i), we have different 2D slices
 *   - a 2D slice (x,y) is represented with the (j,k) dimensions.
 * That is, given an image with 100 slices of size 123 x 456,
 *   ALNSB_PTRTOARRAY3D(data, im, 100, 123, 456);
 * Will declare a variable im, which can be indexed im[i][j][k] with i
 * iterating on the various slices, and j,k iterating inside a slice
 * (j for rows, k for columns).
 *
 * This layout provides slow indexing along slices and fast/stride-one
 * access along pixels which are neighbors along the y dimension in
 * the image.
 *
 */



/**
 *
 *
 *
 */
image3D* image3D_alloc (size_t slices, size_t rows, size_t cols,
			int image_type)
{
  if (slices == 0 || rows == 0 || cols == 0)
    {
      fprintf (stderr, "[%s][ERROR] Cannot allocate an image of size %d x %d x %d\n", ALNSB_PROJECT_STR, slices, rows, cols);
      exit (1);
    }
  image3D* ret = alnsb_calloc (sizeof(image3D), 1);
  ret->data = NULL;
  ret->image_type = image_type;
  ret->image3D = ret;
  ret->slices = slices;
  ret->rows = rows;
  ret->cols = cols;
  size_t nb_pix = ret->slices * ret->rows * ret->cols;
  ret->num_pixels = nb_pix;
  switch (image_type)
    {
    case ALNSB_IMAGE_BINARY: ret->pixel_sz = sizeof(ALNSB_IMAGE_TYPE_BIN); break;
    case ALNSB_IMAGE_INTEGER: ret->pixel_sz =sizeof(ALNSB_IMAGE_TYPE_INT); break;
    case ALNSB_IMAGE_REAL: ret->pixel_sz = sizeof(ALNSB_IMAGE_TYPE_REAL); break;
    default:
      {
      fprintf (stderr, "[%s][ERROR] Cannot allocate an image of type %d (unsupported image type)\n", ALNSB_PROJECT_STR, image_type);
      exit (1);
      }
    }
  ret->data = alnsb_calloc (ret->pixel_sz, nb_pix);
  ret->object_size = ret->pixel_sz * nb_pix;
  ret->ref_counter = 0;

  return ret;
}


/**
 *
 *
 *
 */
void image3D_free (image3D* im)
{
  assert(im != NULL);
  if (!im)
    return;
  if (im->ref_counter > 0)
    im->ref_counter--;
  else
    {
      free (im->data);
      free (im);
    }
}



/**
 *
 *
 *
 */
image3D* image3D_duplicate (image3D* im)
{
  assert(im != NULL);
  image3D* ret = image3D_alloc (im->slices, im->rows, im->cols, im->image_type);
  memcpy (ret->data, im->data, im->object_size);

  return ret;
}



/**
 *
 *
 *
 *
 */
void image3D_copy (image3D* imDst, image3D* imSrc)
{
  if (imSrc == NULL || imDst == NULL)
    {
      fprintf (stderr, "[%s][ERROR] Copying image %p to %p impossible (NULL pointer)\n", ALNSB_PROJECT_STR, imSrc, imDst);
      exit (1);
    }
  if (imSrc->slices != imDst->slices ||
      imSrc->rows != imDst->rows ||
      imSrc->cols != imDst->cols ||
      imSrc->image_type != imDst->image_type)
    {
      fprintf (stderr, "[%s][ERROR] Only copying images of the same size and same type is supported\n", ALNSB_PROJECT_STR);
      exit (1);
    }
  memcpy (imDst->data, imSrc->data, imDst->object_size);
}



/**
 *
 *
 */
extern
image3DBin* image3DBin_alloc(size_t slices, size_t rows, size_t cols)
{
  return (image3DBin*)image3D_alloc (slices, rows, cols, ALNSB_IMAGE_BINARY);
}


/**
 *
 *
 */
image3DInt* image3DInt_alloc(size_t slices, size_t rows, size_t cols)
{
  return (image3DInt*)image3D_alloc (slices, rows, cols, ALNSB_IMAGE_INTEGER);
}

/**
 *
 *
 */
image3DReal* image3DReal_alloc(size_t slices, size_t rows, size_t cols)
{
  return (image3DReal*)image3D_alloc (slices, rows, cols, ALNSB_IMAGE_REAL);
}



/**
 *
 *
 */
void* alnsb_alloc3d (size_t data_sz, size_t num_elt_1, size_t num_elt_2,
		     size_t num_elt_3)
{
  size_t nelt = num_elt_1 * num_elt_2 * num_elt_3;

  return alnsb_calloc (data_sz, nelt);
}


/**
 * Generic image reader. It reads data from a file
 * 'filename' (e.g., "foo.dat"). Arguments are:
 *
 * 'image_type' is one of ALNSB_IMAGE_RAW or ALNSB_IMAGE_TXT. 'RAW'
 * images are data files containing for each pixel the 32/64 bit
 * representation of the pixel value (depending on argument
 * 'data_type). 'TXT' is a plain text file using ASCII to represent
 * numbers, each delimited by space, tab or newline.
 *
 * 'data_type' is one of ALNSB_PIXEL_DATA_TYPE_INT or
 * ALNSB_PIXEL_DATA_TYPE_FLOAT or ALNSB_PIXEL_DATA_TYPE_DOUBLE. It is
 * the type of data used for the pixel representation in the file/in
 * the program.
 *
 * 'nb_elt' is the number of pixels to read.
 *
 */
void* alnsb_readImg (char* filename, int image_type, int data_type,
		     size_t nb_elt)
{
  size_t data_sz = 0;
  switch (data_type)
    {
    case ALNSB_IMAGE_INTEGER: data_sz = sizeof(ALNSB_IMAGE_TYPE_INT); break;
    case ALNSB_IMAGE_REAL: data_sz = sizeof(ALNSB_IMAGE_TYPE_REAL); break;
    case ALNSB_IMAGE_BINARY: data_sz = sizeof(ALNSB_IMAGE_TYPE_BIN); break;
    default: break;
    }
  if (data_sz == 0 || nb_elt == 0)
    {
      fprintf (stderr, "[%s][ERROR] Cannot load image %s of %d elements of size %d bytes\n", ALNSB_PROJECT_STR, filename, data_sz, nb_elt);
      exit (1);
    }

  FILE* f = fopen (filename, "r");
  if (f == NULL)
    {
      fprintf (stderr, "[%s][ERROR] Cannot open file %s\n", ALNSB_PROJECT_STR,
	       filename);
      exit (1);
    }

  void* ret = alnsb_calloc (data_sz, nb_elt);

#define ALNSB_RD_BUFF_SZ 1024
  char buffer[ALNSB_RD_BUFF_SZ];
  size_t nb_read = 0;
  size_t i;
  size_t nc;

  if (data_type == ALNSB_IMAGE_RAW)
    {
      char* cur = ret;
      while ((nc = fread (buffer, sizeof(char), ALNSB_RD_BUFF_SZ, f)))
	{
	  for (i = 0; i < nc; ++i)
	    {
	      if (nb_read < data_sz * nb_elt)
		*(cur++) = buffer[i];
	      ++nb_read;
	    }
	}
      fclose (f);

      if (! feof (f) || nb_read != data_sz * nb_elt)
	{
	  fprintf (stderr, "[%s][ERROR] Read %d bytes of data in file %s (expected %d)%s\n",
		   ALNSB_PROJECT_STR, nb_read, filename, data_sz * nb_elt,
		   feof (f) ? " and EOF reached" : "");
	  exit (1);
	}
    }
  else if (data_type == ALNSB_IMAGE_TXT)
    {
      int* data_i = (int*)ret;
      float* data_f = (float*)ret;
      double* data_d = (double*)ret;
      char b[2];
      for (i = 0; i < nb_elt; ++i)
	{
	  int pos = 0;
	  while ((nc = fread (b, sizeof(char), 1, f)))
	    {
	      if (b[0] != ' ' && b[0] != '\0' && b[0] != '\n' && b[0] != '\t')
		buffer[pos++] = b[0];
	      else
		break;
	    }
	  if (nc == 0)
	    break;
	  if (pos == 0)
	    {
	      // Skip separators.
	      --i;
	      continue;
	    }
	  buffer[pos] = '\0';
	  char* endptr = NULL;
	  int has_error = 0;
	  switch (data_type)
	    {
	    case ALNSB_PIXEL_DATA_TYPE_INT:
	      {
		int val = strtol (buffer, &endptr, 10);
		if (endptr)
		  has_error = 1;
		else
		  *(data_i++) = val;
		break;
	      }
	    case ALNSB_PIXEL_DATA_TYPE_FLOAT:
	      {
		float val = strtof (buffer, &endptr);
		if (endptr)
		  has_error = 1;
		else
		  *(data_f++) = val;
		break;
	      }
	    case ALNSB_PIXEL_DATA_TYPE_DOUBLE:
	      {
		double val = strtod (buffer, &endptr);
		if (endptr)
		  has_error = 1;
		else
		  *(data_d++) = val;
		break;
	      }
	    }
	  if (has_error)
	    {
	      fprintf (stderr, "[%s][ERROR] Cannot convert number \"%s\" at position %d in file %s\n", ALNSB_PROJECT_STR, buffer, i, filename);
	      fclose (f);
	      exit (1);
	    }
	}

      size_t remainder = 0;
      while ((nc = fread (buffer, sizeof(char), ALNSB_RD_BUFF_SZ, f)))
	remainder += nc;
      fclose (f);
      if (remainder > 0)
	{
	  fprintf (stderr, "[%s][ERROR] File %s has %d pixels (expected %d)\n", ALNSB_PROJECT_STR, filename, i + remainder, nb_elt);
	  exit (1);
	}
    }
  else
    {
      fprintf (stderr, "[%s][ERROR] Unsupported image type: %d\n",
	       ALNSB_PROJECT_STR, data_type);
      fclose (f);
      exit (1);
    }

  return ret;
}


/**
 *
 *
 *
 */
void alnsb_writeImg (void* data, char* filename, int image_type, int data_type,
		     size_t nb_elt)
{


}



/**
 *
 *
 *
 */
void alnsb_displayImg (void* data, char* viewscript, int image_type, int data_type,
		     size_t nb_elt)
{


}
