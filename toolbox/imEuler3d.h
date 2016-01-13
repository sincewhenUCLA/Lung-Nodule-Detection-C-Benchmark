/**
 * imEuler3d.h: this file is part of the ALNSB project.
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

#ifndef ALNSB_TOOLBOX_IMEULER3D_H
# define ALNSB_TOOLBOX_IMEULER3D_H

# include <utilities/images.h>

extern
int
alnsb_imEuler3d_bin3d(ALNSB_IMAGE_TYPE_BIN* __ALNSB_RESTRICT_PTR im3d,
		      int dim0, int dim1, int dim2);


#endif // !ALNSB_TOOLBOX_IMEULER3D_H
