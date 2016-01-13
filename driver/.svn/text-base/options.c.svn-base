/**
 * options.c: this file is part of the ALNSB project.
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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <strings.h>

#include <utilities/common.h>
#include <utilities/environment.h>

struct s_options
{
  char* longopt;
  char* shortopt;
  int has_arg;
  void* dst_ptr;
  int dst_ptr_type;
  char* description;
};


#define LSCAD_OPT_NONE 0
#define LSCAD_OPT_INT  1
#define LSCAD_OPT_REAL 2
#define LSCAD_OPT_STR  3

void read_value (char* arg, void* dst, int type)
{
  if (dst == NULL)
    return;
  switch (type)
    {
    case LSCAD_OPT_NONE:
      {
	int* dstptr = (int*) dst;
	*dstptr = 1;
	break;
      }
    case LSCAD_OPT_INT:
      {
	int* dstptr = (int*) dst;
	int val = 1;
	if (arg != NULL)
	  val = atoi (arg);
	*dstptr = val;
	break;
      }
    case LSCAD_OPT_REAL:
      {
	double* dstptr = (double*) dst;
	double val = 1;
	if (arg != NULL)
	  val = atof (arg);
	*dstptr = val;
	break;
      }
    case LSCAD_OPT_STR:
      {
	char** dstptr = (char**) dst;
	if (arg != NULL)
	  *dstptr = strdup (arg);
	break;
      }
    default: break;
    }
}

void print_value (FILE* f, void* src, int type)
{
  if (src == NULL)
    return;
  switch (type)
    {
    case LSCAD_OPT_NONE: 
    case LSCAD_OPT_INT:
      {
	int* srcptr = (int*) src;
	fprintf (f, "%d", *srcptr);
	break;
      }
    case LSCAD_OPT_REAL:
      {
	double* srcptr = (double*) src;
	fprintf (f, "%f", *srcptr);
	break;
      }
    case LSCAD_OPT_STR:
      {
	char** srcptr = (char**) src;
	fprintf (f, "%s", *srcptr);
	break;
      }
    default: break;
    }
}

int alnsb_getopts (int argc, char** argv, s_alnsb_environment_t* env)
{
  // Declare all options.
  struct s_options myopts[] = {
    // Help must always be first option.
    { "--help", "-h", 0, NULL, LSCAD_OPT_NONE },

    // Others can be in any order, but help will print them in this order.
    { "--enable-timing", NULL, 0, &(env->timer), LSCAD_OPT_NONE,
      "[general] Time each pass" },
    { "--verbose-level", NULL, 1, &(env->verbose_level), LSCAD_OPT_INT,
      "[general] verbosity level" },
    { "--show-environment", NULL, 0, &(env->show_environment), LSCAD_OPT_NONE,
      "[general] Show all program parameter values" },

    { "--patient_folder", "-pf", 1, &(env->patient_name), LSCAD_OPT_STR,
      "[images] Name of folder in images/ containing the images" },

    { "--num_slices", "-nz", 1, &(env->num_slices), LSCAD_OPT_INT,
      "[images] Number of 2D slices" },
    { "--slice-size-x", "-nx", 1, &(env->slice_size_x), LSCAD_OPT_INT,
      "[images] Number of pixels in the 'x' dimension of a 2D slice" },
    { "--slice-size-y", "-ny", 1, &(env->slice_size_y), LSCAD_OPT_INT,
      "[images] Number of pixels in the 'y' dimension of a 2D slice" },
    { "--rotate-90-right", NULL, 0, &(env->transpose_input), LSCAD_OPT_NONE,
      "[images] Rotate 90 degrees right" },

    { "--scanner_pixel_spacing_x_mm", NULL, 1, &(env->scanner_pixel_spacing_x_mm), LSCAD_OPT_REAL,
      "[scanner] Pixel spacing in 'x' dimension, in mm" },
    { "--scanner_pixel_spacing_y_mm", NULL, 1, &(env->scanner_pixel_spacing_y_mm), LSCAD_OPT_REAL,
      "[scanner] Pixel spacing in 'y' dimension, in mm" },
    { "--scanner_slice_thickness_mm", NULL, 1, &(env->scanner_slice_thickness_mm), LSCAD_OPT_REAL,
      "[scanner] Slice thickness, in mm" },

    { "--thresold_upperBand", NULL, 1, &(env->thresold_upperBand), LSCAD_OPT_INT,
      "[thresolding] Maximal pixel value" },
    { "--thresold_lowerBand", NULL, 1, &(env->thresold_lowerBand), LSCAD_OPT_INT,
      "[thresolding] Minimal pixel value" },

    { "--thresolding_skip4", NULL, 0, &(env->segmentationMask_skip_4_slices), LSCAD_OPT_NONE,
      "[thresolding] Skip first/last 4 slices" },

    { "--segmentation_lp", NULL, 1, &(env->segmentation_lp), LSCAD_OPT_REAL,
      "[segmentation] lp parameter value" },
    { "--segmentation_errb0", NULL, 1, &(env->segmentation_errb[0]), LSCAD_OPT_REAL,
      "[segmentation] errb[0] parameter value" },
    { "--segmentation_errb1", NULL, 1, &(env->segmentation_errb[1]), LSCAD_OPT_REAL,
      "[segmentation] errb[1] parameter value" },
    { "--segmentation_ulab0", NULL, 1, &(env->segmentation_ulab[0]), LSCAD_OPT_REAL,
      "[segmentation] ulab[0] parameter value" },
    { "--segmentation_ulab1", NULL, 1, &(env->segmentation_ulab[1]), LSCAD_OPT_REAL,
      "[segmentation] ulab[1] parameter value" },

    { "--segmentation_cc", NULL, 1, &(env->segmentation_cc), LSCAD_OPT_REAL,
      "[segmentation] cc parameter value" },
    { "--segmentation_c_convergence", NULL, 1, &(env->segmentation_c_convergence), LSCAD_OPT_REAL,
      "[segmentation] c_convergence parameter value" },
    { "--segmentation_steps", NULL, 1, &(env->segmentation_steps), LSCAD_OPT_REAL,
      "[segmentation] steps parameter value" },
    { "--segmentation_beta", NULL, 1, &(env->segmentation_beta), LSCAD_OPT_REAL,
      "[segmentation] beta parameter value" },
    { "--segmentation_err_c", NULL, 1, &(env->segmentation_err_c), LSCAD_OPT_REAL,
      "[segmentation] err_c parameter value" },
    { "--segmentation_max_err", NULL, 1, &(env->segmentation_max_err), LSCAD_OPT_REAL,
      "[segmentation] max_err parameter value" },
    { "--segmentation_max_steps", NULL, 1, &(env->segmentation_max_steps), LSCAD_OPT_INT,
      "[segmentation] max_steps parameter value" },

    { "--preselection_diameterMin", NULL, 1, &(env->preselection_diameterMin), LSCAD_OPT_REAL,
      "[preselection] Minimal nodule diameter" },
    { "--preselection_diameterMax", NULL, 1, &(env->preselection_diameterMax), LSCAD_OPT_REAL,
      "[preselection] Maximal nodule diameter" },
    { "--preselection_elongationMax", NULL, 1, &(env->preselection_elongationMax), LSCAD_OPT_REAL,
      "[preselection] Maximal nodule elongation" },
    { "--preselection_circulMin", NULL, 1, &(env->preselection_circulMin), LSCAD_OPT_REAL,
      "[preselection] Minimal circularity" },

    { "--classifier_num_features", NULL, 1, &(env->classifier_num_features), LSCAD_OPT_INT,
      "[classifier] Number of features" },
    { "--classifier_posMat_filename", NULL, 1, &(env->classifier_positive_featMat_filename), LSCAD_OPT_STR,
      "[classifier] Name of file containing positive samples" },
    { "--classifier_negMat_filename", NULL, 1, &(env->classifier_negative_featMat_filename), LSCAD_OPT_STR,
      "[classifier] Name of file containing negative samples" },
    { "--classifier_stdVec_filename", NULL, 1, &(env->classifier_stdFeat_filename), LSCAD_OPT_STR,
      "[classifier] Name of file containing standard feature vector" },
    { "--classifier_meanVec_filename", NULL, 1, &(env->classifier_meanFeat_filename), LSCAD_OPT_STR,
      "[classifier] Name of file containing standard feature vector" },
    { "--classifier_active_features", NULL, 1, &(env->classifier_featvect_str), LSCAD_OPT_STR,
      "[classifier] Feature vector (e.g., {0,1,0} with no space)" },

    { "--emtv-skip", NULL, 0, &(env->pass_options[0].load_pass_result), LSCAD_OPT_NONE,
      "[emtv] Load saved pass result instead of executing it (inactive)" },

    { "--emtv-display", NULL, 0, &(env->pass_options[0].display_pass_result), LSCAD_OPT_NONE,
      "[emtv] Display pass result" },

    { "--rotation-skip", NULL, 0, &(env->pass_options[1].load_pass_result), LSCAD_OPT_NONE,
      "[rotation] Load saved pass result instead of executing it" },
    { "--rotation-display", NULL, 0, &(env->pass_options[1].display_pass_result), LSCAD_OPT_NONE,
      "[rotation] Display pass result" },

    { "--levelscale-skip", NULL, 0, &(env->pass_options[2].load_pass_result), LSCAD_OPT_NONE,
      "[levelscale] Load saved pass result instead of executing it" },
    { "--levelscale-display", NULL, 0, &(env->pass_options[2].display_pass_result), LSCAD_OPT_NONE,
      "[levelscale] Display pass result" },

    { "--segmentation-skip", NULL, 0, &(env->pass_options[3].load_pass_result), LSCAD_OPT_NONE,
      "[segmentation] Load saved pass result instead of executing it" },
    { "--segmentation-display", NULL, 0, &(env->pass_options[3].display_pass_result), LSCAD_OPT_NONE,
      "[segmentation] Display pass result" },
    { "--segmentationMask-skip", NULL, 0, &(env->pass_options[4].load_pass_result), LSCAD_OPT_NONE,
      "[segmentationMask] Load saved pass result instead of executing it" },
    { "--segmentationMask-display", NULL, 0, &(env->pass_options[4].display_pass_result), LSCAD_OPT_NONE,
      "[segmentationMask] Display pass result" },
    { "--preselection-skip", NULL, 0, &(env->pass_options[5].load_pass_result), LSCAD_OPT_NONE,
      "[preselection] Load saved pass result instead of executing it" },
    { "--preselection-display", NULL, 0, &(env->pass_options[5].display_pass_result), LSCAD_OPT_NONE,
      "[preselection] Display pass result" },
    { "--featureExtraction-skip", NULL, 0, &(env->pass_options[6].load_pass_result), LSCAD_OPT_NONE,
      "[featureExtraction] Load saved pass result instead of executing it" },
    { "--featureExtraction-display", NULL, 0, &(env->pass_options[6].display_pass_result), LSCAD_OPT_NONE,
      "[featureExtraction] Display pass result" },
    { "--classification-skip", NULL, 0, &(env->pass_options[7].load_pass_result), LSCAD_OPT_NONE,
      "[classification] Load saved pass result instead of executing it" },
    { "--classification-display", NULL, 0, &(env->pass_options[7].display_pass_result), LSCAD_OPT_NONE,
      "[classification] Display pass result" },
    
    // Terminator. Must be null for longopt.
    { NULL, NULL, 0, NULL, LSCAD_OPT_NONE}
  };

  assert(env != NULL);
  
  int i, j;
  for (i = 1; i < argc; ++i)
    {
      for (j = 0; myopts[j].longopt != NULL; ++j)
	{
	  if (! strcmp (argv[i], myopts[j].longopt) ||
	      (myopts[j].shortopt && ! strcmp (argv[i], myopts[j].shortopt)))
	    {
	      if (j == 0)
		{
		  fprintf (stderr, "[ALNSB] List of available options:\n");
		  int k;
		  for (k = 1; myopts[k].longopt != NULL; ++k)
		    {
		      int sk = 4 - strlen (myopts[k].longopt) / 8;
		      int l;
		      fprintf (stderr, "%s ", myopts[k].longopt);
		      fprintf (stderr, "%s",
			       (myopts[k].has_arg ? "<value-needed> " : "               "));
		      for (l = 0; l < sk; ++l)
			fprintf (stderr, "\t");
		      fprintf (stderr, "%s\n", myopts[k].description);
		    }
		  exit (1);
		}
	      if (! myopts[j].has_arg)
		{
		  read_value (NULL, myopts[j].dst_ptr, myopts[j].dst_ptr_type);
		  break;
		}
	      else if (myopts[j].has_arg && i < argc - 1)
		{
		  read_value (argv[i + 1], myopts[j].dst_ptr, 
			      myopts[j].dst_ptr_type);
		  ++i;
		  break;
		}
	      else
		fprintf (stderr, "[ALNSB][ERROR] Missing argument to %s\n",
			 argv[i]);
	    }
	}
      if (myopts[j].longopt == NULL)
	{
	  fprintf (stderr, "[ALNSB][ERROR] Unsupported option %s %d %d\n",
		   argv[i], i, argc);
	  exit (1);
	}
    }
  if (env->verbose_level > 1 || env->show_environment)
    {
      fprintf (stderr, "[ALNSB][DEBUG] Environment values:\n");
      for (i = 1; myopts[i].longopt != NULL; ++i)
	{
	  fprintf (stderr, "%s", myopts[i].longopt);
	  int sk = 4 - strlen (myopts[i].longopt) / 8;
	  for (j = 0; j < sk; ++j)
	    fprintf (stderr, "\t");
	  fprintf (stderr, "=>\t");
	  print_value (stderr, myopts[i].dst_ptr, myopts[i].dst_ptr_type);
	  fprintf (stderr, "\n");
	}
    }
}
