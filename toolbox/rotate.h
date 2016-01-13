/**
 * rotate.h: this file is part of the ALNSB project.
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

#ifndef ALNSB_TOOLBOX_ROTATE_H
# define ALNSB_TOOLBOX_ROTATE_H

# include <utilities/images.h>

extern
void alnsb_rotate_slices_real_img (image3DReal* __ALNSB_RESTRICT_PTR img);



#endif // !ALNSB_TOOLBOX_ROTATE_H
