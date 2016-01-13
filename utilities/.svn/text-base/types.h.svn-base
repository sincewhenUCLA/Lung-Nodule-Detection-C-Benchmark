/**
 * types.h: this file is part of the ALNSB project.
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
#ifndef ALNSB_TYPES_H
# define ALNSB_TYPES_H


/**
 * Image type. Set to TXT if the image is a text file containing
 * numbers, in slice-major-then-row-major order.  Set to RAW if the
 * image is a raw data file of 64-bit double-precision numbers, in the
 * slice-major-then-row-major order.
 */
#define ALNSB_IMAGE_TXT 1
#define ALNSB_IMAGE_RAW 2


/**
 * Image element scalar type.
 *
 */
#define ALNSB_PIXEL_DATA_TYPE_INT	1
#define ALNSB_PIXEL_DATA_TYPE_FLOAT	2
#define ALNSB_PIXEL_DATA_TYPE_DOUBLE	3

/**
 * Never more than 10 stages (currently: 6)
 *
 */
#define ALNSB_MAX_NUMBER_OF_PIPELINE_STAGES 10


/**
 * Image kind.
 *
 */
#define ALNSB_IMAGE_BINARY	1  // will be supported in the future. See them as integer images for the moment.
#define ALNSB_IMAGE_INTEGER	2
#define ALNSB_IMAGE_REAL	3


/**
 * Type of data for the feature matrices. Default to double.
 */
#define ALNSB_MAT_TYPE		double
/**
 * Type of data for the real images. Default to double.
 */
#define ALNSB_IMAGE_TYPE_REAL	float
/**
 * Type of data for the binary images. Default to int.
 */
#define ALNSB_IMAGE_TYPE_BIN	float
/**
 * Type of data for the integer images. Default to int.
 */
#define ALNSB_IMAGE_TYPE_INT	int
/**
 * Type of data for the txt images. Default to char.
 */
#define ALNSB_IMAGE_TYPE_TXT	char



#endif // !ALNSB_TYPES_H
