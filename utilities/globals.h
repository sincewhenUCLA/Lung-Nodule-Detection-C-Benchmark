/**
 * globals.h: this file is part of the ALNSB project.
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
#ifndef ALNSB_GLOBALS_H
# define ALNSB_GLOBALS_H

/**
 * Project name, used for console printing.
 */
#define ALNSB_PROJECT_STR "ALNSB"

/**
 * Image type. Set to TXT if the image is a text file containing
 * numbers, in slice-major-then-row-major order.  Set to RAW if the
 * image is a raw data file of 64-bit double-precision numbers, in the
 * slice-major-then-row-major order.
 */
#define ALNSB_IMAGE_TXT 1
#define ALNSB_IMAGE_RAW 2


/**
 * Type of data for the feature matrices. Default to double.
 */
#define ALNSB_MAT_TYPE		double
/**
 * Type of data for the real images. Default to double.
 */
#define ALNSB_IMAGE_TYPE_REAL	double
/**
 * Type of data for the binary images. Default to int.
 */
#define ALNSB_IMAGE_TYPE_BIN	int



#endif // !ALNSB_TYPES_H
