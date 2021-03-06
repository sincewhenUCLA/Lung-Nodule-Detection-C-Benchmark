/**
 * environment.c: this file is part of the ALNSB project.
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
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>

#include <utilities/environment.h>
#include <utilities/memfuncs.h>
#include <utilities/images.h>


s_alnsb_environment_t* alnsb_environment_malloc ()
{
  s_alnsb_environment_t* env = alnsb_calloc (sizeof(s_alnsb_environment_t), 1);

  // Initialize default values.

  // Input image info.
  env->num_slices = 124;
  env->slice_size_x = 716;
  env->slice_size_y = 716;

  // Scanner info.
  env->scanner_pixel_spacing_x_mm = 0.7;
  env->scanner_pixel_spacing_y_mm = 0.7;
  env->scanner_slice_thickness_mm = 1.25;

  env->thresold_upperBand = 80;
  env->thresold_lowerBand = 0;
  env->segmentation_lp = 1e-13;
  env->segmentation_errb[0] = 1e-1;
  env->segmentation_errb[1] = 5e-4;
  env->segmentation_ulab[0] = 0.1;
  env->segmentation_ulab[1] = 0.4;

  env->segmentation_cc = 0.35;
  env->segmentation_c_convergence = 3e-4;
  env->segmentation_steps = 0.11;
  env->segmentation_beta = 2;
  env->segmentation_err_c = 1;
  env->segmentation_max_err = 0.0;
  env->segmentation_max_steps = 300;

  env->preselection_diameterMin = 5;
  env->preselection_diameterMax = 30;
  env->preselection_elongationMax = 4;
  env->preselection_circulMin = 1.0/6.0;

  env->verbose_level = 0;
  env->show_environment = 0;

  // Classifier info.
  env->classifier_num_features = 27;
  env->classifier_positive_featMat_filename = strdup ("data/classifier-1/SelectedPositiveSamples.dat");
  env->classifier_positive_featMat_numSamples = 1;
  env->classifier_negative_featMat_filename = strdup ("data/classifier-1/SelectedNegativeSamples.dat");
  env->classifier_negative_featMat_numSamples = 21;
  env->classifier_stdFeat_filename = strdup ("data/classifier-1/stdFeature.dat");
  env->classifier_meanFeat_filename = strdup ("data/classifier-1/meanFeature.dat");

  env->patient_name = strdup("NLST_R0960B_OUT4");

  env->dump_images = 1;

  env->pass_options[0].pass_name = "emtv";
  env->pass_options[0].load_pass_result = 1;
  env->pass_options[0].display_pass_result = 0;
  env->pass_options[1].pass_name = "rotation";
  env->pass_options[1].load_pass_result = 0;
  env->pass_options[1].display_pass_result = 0;
  env->pass_options[2].pass_name = "levelscale";
  env->pass_options[2].load_pass_result = 0;
  env->pass_options[2].display_pass_result = 0;
  env->pass_options[3].pass_name = "segmentation";
  env->pass_options[3].load_pass_result = 0;
  env->pass_options[3].display_pass_result = 0;
  env->pass_options[4].pass_name = "segmentationMask";
  env->pass_options[4].load_pass_result = 0;
  env->pass_options[4].display_pass_result = 0;
  env->pass_options[5].pass_name = "preselection";
  env->pass_options[5].load_pass_result = 0;
  env->pass_options[5].display_pass_result = 0;
  env->pass_options[6].pass_name = "featureExtraction";
  env->pass_options[6].load_pass_result = 0;
  env->pass_options[6].display_pass_result = 0;
  env->pass_options[7].pass_name = "classification";
  env->pass_options[7].load_pass_result = 0;
  env->pass_options[7].display_pass_result = 0;


  env->classifier_featvect_str =
    //strdup ("{0,1,1,1,0,1,1,0,1,0,1,1,0,1,0,0,0,0,0,0,1,1,1,0,1,1,0}");
  strdup ("{1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}");
  env->classifier_active_features =
    (char*) malloc (sizeof(char) * env->classifier_num_features);
  char* feats = env->classifier_featvect_str;
  ++feats;
  int i;
  for (i = 0; i < env->classifier_num_features; ++i)
    {
      env->classifier_active_features[i] = (*feats == '0' ? 0 : 1);
      feats += 2;
    }

  return env;
}


void alnsb_environment_print (s_alnsb_environment_t* env)
{
  if (env == NULL)
    {
      fprintf (stderr, "[INFO] Environment is null\n");
      return;
    }

}

void alnsb_environment_free (s_alnsb_environment_t* env)
{
  if (env)
    free (env);

}
