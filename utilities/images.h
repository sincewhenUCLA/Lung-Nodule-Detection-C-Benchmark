/**
 * images.h: this file is part of the ALNSB project.
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
#ifndef ALNSB_IMAGES_H
# define ALNSB_IMAGES_H

#include <stdio.h>
#include <stdlib.h>
#include <utilities/common.h>
#include <utilities/environment.h>
#include <utilities/types.h>

union ls_image_type {
  ALNSB_IMAGE_TYPE_BIN* img_bin;
  ALNSB_IMAGE_TYPE_INT* img_int;
  ALNSB_IMAGE_TYPE_REAL* img_real;
  ALNSB_IMAGE_TYPE_TXT* img_txt;
};
typedef union ls_image_type ls_image_type_t;

struct s_image3D
{
  // Pointer to self class type (convenience)
  struct s_image3D* image3D;
  // Number of pixel in each dimension.
  size_t slices, rows, cols;
  // Total number of pixels, = slices * rows * cols (convenience).
  size_t num_pixels;
  // Image type can be real or int or txt.
  // ALNSB_IMAGE_BINARY | ALNSB_IMAGE_INT | ALNSB_IMAGE_REAL | ALNSB_IMAGE_TXT
  int image_type;
  // size of a pixel, in byte.
  size_t pixel_sz;
  // Number of bytes of the matrix, may differ from
  // num_pixels*pixel_sz for txt images.
  size_t object_size;
  // Data.
  void* data;
  // Reference counter.
  int ref_counter;
};
typedef struct s_image3D image3D;

struct s_image3DReal
{
  // Pointer to super class type (convenience)
  image3D* image3D;
  // Number of pixel in each dimension.
  size_t slices, rows, cols;
  // Total number of pixels, = slices * rows * cols (convenience).
  size_t num_pixels;
  // Image type can be real or int or txt.
  // ALNSB_IMAGE_BINARY | ALNSB_IMAGE_INT | ALNSB_IMAGE_REAL | ALNSB_IMAGE_TXT
  int image_type;
  // size of a pixel, in byte.
  size_t pixel_sz;
  // Number of bytes of the matrix, may differ from
  // num_pixels*pixel_sz for txt images.
  size_t object_size;
  // Typed data.
  ALNSB_IMAGE_TYPE_REAL* data;
  // Reference counter.
  int ref_counter;
};
typedef struct s_image3DReal image3DReal;


struct s_image3DInt
{
  // Pointer to super class type (convenience)
  image3D* image3D;
  // Number of pixel in each dimension.
  size_t slices, rows, cols;
  // Total number of pixels, = slices * rows * cols (convenience).
  size_t num_pixels;
  // Image type can be real or int or txt.
  // ALNSB_IMAGE_BINARY | ALNSB_IMAGE_INT | ALNSB_IMAGE_REAL | ALNSB_IMAGE_TXT
  int image_type;
  // size of a pixel, in byte.
  size_t pixel_sz;
  // Number of bytes of the matrix, may differ from
  // num_pixels*pixel_sz for txt images.
  size_t object_size;
  // Typed data.
  ALNSB_IMAGE_TYPE_INT* data;
  // Reference counter.
  int ref_counter;
};
typedef struct s_image3DInt image3DInt;



struct s_image3DBin
{
  // Pointer to super class type (convenience)
  image3D* image3D;
  // Number of pixel in each dimension.
  size_t slices, rows, cols;
  // Total number of pixels, = slices * rows * cols (convenience).
  size_t num_pixels;
  // Image type can be real or int or txt.
  // ALNSB_IMAGE_BINARY | ALNSB_IMAGE_INT | ALNSB_IMAGE_REAL | ALNSB_IMAGE_TXT
  int image_type;
  // size of a pixel, in byte.
  size_t pixel_sz;
  // Number of bytes of the matrix, may differ from
  // num_pixels*pixel_sz for txt images.
  size_t object_size;
  // Typed data.
  ALNSB_IMAGE_TYPE_BIN* data;
  // Reference counter.
  int ref_counter;
};
typedef struct s_image3DBin image3DBin;


/**
 *
 *
 *
 */
extern
image3D* image3D_alloc (size_t slices, size_t rows, size_t cols,
			int image_type);

/**
 *
 *
 */
extern
image3DBin* image3DBin_alloc(size_t slices, size_t rows, size_t cols);

/**
 *
 *
 */
extern
image3DInt* image3DInt_alloc(size_t slices, size_t rows, size_t cols);

/**
 *
 *
 */
extern
image3DReal* image3DReal_alloc(size_t slices, size_t rows, size_t cols);

/**
 *
 *
 *
 */
extern
void image3D_free (image3D* im);


/**
 *
 *
 *
 */
extern
image3D* image3D_duplicate (image3D* im);



/**
 *
 *
 *
 *
 */
extern
void image3D_copy (image3D* imDst, image3D* imSrc);



#define ALNSB_IMBinTo3D(imagePtr,arrname)				\
  ALNSB_IMAGE_TYPE_BIN (*arrname)[(imagePtr)->rows][(imagePtr)->cols] = (ALNSB_IMAGE_TYPE_BIN (*)[(imagePtr)->rows][(imagePtr)->cols])((imagePtr)->data)
#define ALNSB_IMIntTo3D(imagePtr,arrname)				\
  ALNSB_IMAGE_TYPE_INT (*arrname)[(imagePtr)->rows][(imagePtr)->cols] = (ALNSB_IMAGE_TYPE_INT (*)[(imagePtr)->rows][(imagePtr)->cols])((imagePtr)->data)

#define ALNSB_IMRealTo3D(imagePtr,arrname)				\
  ALNSB_IMAGE_TYPE_REAL (*arrname)[(imagePtr)->rows][(imagePtr)->cols] = (ALNSB_IMAGE_TYPE_REAL (*)[(imagePtr)->rows][(imagePtr)->cols])((imagePtr)->data)

#define ALNSB_IMRealTo1D(imagePtr,arrname)				\
  ALNSB_IMAGE_TYPE_REAL* __ALNSB_RESTRICT_PTR arrname = (imagePtr)->data

#define ALNSB_IMBinTo1D(imagePtr,arrname)				\
  ALNSB_IMAGE_TYPE_BIN* __ALNSB_RESTRICT_PTR arrname = (imagePtr)->data

#define ALNSB_IMIntTo1D(imagePtr,arrname)				\
  ALNSB_IMAGE_TYPE_INT* __ALNSB_RESTRICT_PTR arrname = (imagePtr)->data


/**
 *
 * Data layout: a 3D image is allocated as a 1-dimensional vector of
 * data, using the macro ALNSB_ALLOC3D(type, size1, size2, size3)
 * which returns a type* pointer.
 *
 * This pointer can be viewed as a 3D image im[i][j][k] using the
 * macro ALNSB_PTRTOARRAY3D(ptrname, arrname, type, size1, size2, size3)
 * to declare a new variable arrname which is a C99 3D array (e.g.,
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



#define ALNSB_ALLOC3D(type, size1, size2, size3) \
  alnsb_alloc3d (sizeof(type), (size1), (size2), (size3))

#define ALNSB_PTRTOARRAY3D(ptrname, arrname, type, size1, size2, size3) \
  type (*arrname)[(size2)][(size3)] = (type (*)[(size2)][(size3)])(ptrname)


#define ALNSB_REAL3D_alloc(im,slices,rows,cols)				\
  image3DReal* __ALNSB_RESTRICT_PTR im##_im = image3DReal_alloc (slices, rows, cols); \
  ALNSB_IMRealTo3D(im##_im,im);

#define ALNSB_REAL3D_free(im)				\
  image3D_free (im##_im->image3D);



#define ALNSB_BIN1D_alloc(im,slices,rows,cols)				\
  image3DBin* __ALNSB_RESTRICT_PTR im##_im = image3DBin_alloc (slices, rows, cols); \
  ALNSB_IMBinTo1D(im##_im,im);


#define ALNSB_REAL1D_alloc(im,slices,rows,cols)				\
  image3DReal* __ALNSB_RESTRICT_PTR im##_im = image3DReal_alloc (slices, rows, cols); \
  ALNSB_IMRealTo1D(im##_im,im);



#define ALNSB_BIN1D_free(im)				\
  image3D_free (im##_im->image3D);

#define ALNSB_2Dslice_from_Bin1D(slice2d, im1d, rows, cols)	\
      ALNSB_IMAGE_TYPE_BIN (*slice2d)[cols] =			\
	(ALNSB_IMAGE_TYPE_BIN (*)[cols])(im1d);



/**
 *
 *
 */
extern
void* alnsb_alloc3d (size_t data_sz, size_t num_elt_1, size_t num_elt_2,
		     size_t num_elt_3);


/**
 *
 *
 */
extern
void* alnsb_readImg (char* filename, int image_type, int data_type,
		     size_t nb_elt);

/**
 *
 *
 *
 */
extern
void alnsb_writeImg (void* data, char* filename, int image_type, int data_type,
		     size_t nb_elt);


#endif // !ALNSB_IMAGES_H
