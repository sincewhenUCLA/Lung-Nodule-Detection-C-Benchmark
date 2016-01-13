/**
 * environment.h: this file is part of the ALNSB project.
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
#ifndef ALNSB_ENVIRONMENT_H
# define ALNSB_ENVIRONMENT_H

# include <utilities/types.h>
# include <utilities/images.h>

// Define the total number of phases in the algorithm (segmentation,
// segmentationMask, preselection, featureExtraction,
// classification). It is used to size the various per-pass
// environment variables such as image input/output filenames, etc.
#define ALNSB_MAX_NUMBER_OF_PHASES 10

struct pass_opts
{
  char* pass_name;
  int   load_pass_result; // set to 1 to not execute the pass and load
			  // its result if available, otherwise will
			  // run the pass.
  int   display_pass_result; // set to 1 to display the pass output.
};

# define EMTV_PASS		0
# define ROTATION_PASS		1
# define LEVELSCALE_PASS	2
# define SEGMENTATION_PASS	3
# define SEGMENTATIONMASK_PASS	4
# define PRESELECTION_PASS	5
# define FEATUREEXTRACTION_PASS 6
# define CLASSIFICATION_PASS	7


struct alnsb_environment {

  // General.
  char*			patient_name;
  int			verbose_level;
  int			show_environment;
  int			timer;

  // Input image info.
  size_t		num_slices;
  size_t		slice_size_x;
  size_t		slice_size_y;
  int			input_image_format[ALNSB_MAX_NUMBER_OF_PHASES]; // ALNSB_IMAGE_TXT | ALNSB_IMAGE_RAW
  char*			input_image_filename[ALNSB_MAX_NUMBER_OF_PHASES];

  // Scanner info.
  double		scanner_pixel_spacing_x_mm;
  double		scanner_pixel_spacing_y_mm;
  double		scanner_slice_thickness_mm;

  // Preparation (rotate==transpose).
  int			transpose_input;
  double		thresold_upperBand;
  double		thresold_lowerBand;

  // Segmentation info.
  double		segmentation_lp;
  double		segmentation_errb[2];
  double		segmentation_ulab[2];
  double		segmentation_cc;
  double		segmentation_c_convergence;
  double		segmentation_steps;
  double		segmentation_beta;
  double		segmentation_err_c;
  double		segmentation_max_err;
  int			segmentation_max_steps;

  // SegmentationMask options.
  int			segmentationMask_skip_4_slices;

  // Preselection options.
  double       		preselection_diameterMin;
  double       		preselection_diameterMax;
  double       		preselection_elongationMax;
  double		preselection_circulMin;

  // Classifier info.
  size_t		classifier_num_features;
  char*			classifier_active_features;
  char*			classifier_positive_featMat_filename;
  size_t		classifier_positive_featMat_numSamples;
  char*			classifier_negative_featMat_filename;
  size_t		classifier_negative_featMat_numSamples;
  char*			classifier_meanFeat_filename;
  char*			classifier_stdFeat_filename;
  char*			classifier_featvect_str;

  // Per-pass info (name, skip/execute, display).
  struct pass_opts	pass_options[ALNSB_MAX_NUMBER_OF_PHASES];

  // Internal use only, must be set to 1.
  int			dump_images;
};
typedef struct alnsb_environment s_alnsb_environment_t;


extern
s_alnsb_environment_t* alnsb_environment_malloc ();

extern
void alnsb_environment_free (s_alnsb_environment_t* env);

extern
void alnsb_print (s_alnsb_environment_t* env);


#endif //!ALNSB_ENVIRONMENT_H
