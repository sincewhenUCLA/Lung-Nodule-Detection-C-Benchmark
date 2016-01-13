/**
 * featureExtraction_step.h: this file is part of the ALNSB project.
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
#ifndef ALNSB_FEATUREEXTRACTION_STEP_H
# define ALNSB_FEATUREEXTRACTION_STEP_H

# include <utilities/types.h>
# include <utilities/images.h>
# include <utilities/environment.h>


extern
void featureExtraction_cpu (s_alnsb_environment_t* __ALNSB_RESTRICT_PTR env,
			    image3DReal* __ALNSB_RESTRICT_PTR inputPrep,
			    image3DBin* __ALNSB_RESTRICT_PTR inputPresel,
			    image3DReal** __ALNSB_RESTRICT_PTR outputFeatures);




#endif //!ALNSB_FEATUREEXTRACTION_STEP_H
