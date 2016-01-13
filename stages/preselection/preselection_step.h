/**
 * preselection_step.h: this file is part of the ALNSB project.
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
#ifndef ALNSB_PRESELECTION_STEP_H
# define ALNSB_PRESELECTION_STEP_H

# include <utilities/types.h>
# include <utilities/images.h>
# include <utilities/environment.h>


extern
void preselection_cpu (s_alnsb_environment_t* __ALNSB_RESTRICT_PTR env,
		       image3DBin* __ALNSB_RESTRICT_PTR input,
		       image3DBin** __ALNSB_RESTRICT_PTR output);




#endif //!ALNSB_PRESELECTION_STEP_H
