/**
 * file_io.h: this file is part of the ALNSB project.
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
#ifndef ALNSB_FILE_IO_H
# define ALNSB_FILE_IO_H
#include <stdio.h>

extern
void* alnsb_read_data_from_binary_file (char* filename, size_t elt_sz,
					size_t elt_count);

extern
void* alnsb_read_data_from_binary_file_nosz (char* filename,
					     size_t elt_sz,
					     size_t max_elt_count,
					     size_t* num_elt_read);

extern
void alnsb_save_data_to_file (void* data, char* filename, size_t elt_sz,
			      size_t elt_count);

extern
void alnsb_save_data_to_file_ascii (float* data, char* filename,
				    size_t elt_count);


#endif // !ALNSB_FILE_IO_H
