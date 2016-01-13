/**
 * pipeline.c: this file is part of the ALNSB project.
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
#include <stdio.h>

#include <driver/pipeline.h>
#include <utilities/memfuncs.h>
#include <utilities/images.h>
#include <utilities/step.h>
#include <utilities/file_io.h>
#include <utilities/timer.h>

#include <stages/rotation/rotation_step.h>
#include <stages/levelscale/levelscale_step.h>
#include <stages/segmentation/segmentation_step.h>
#include <stages/segmentationMask/segmentationMask_step.h>
#include <stages/preselection/preselection_step.h>
#include <stages/featureExtraction/featureExtraction_step.h>
#include <stages/classification/classification_step.h>

static
void display_image (s_alnsb_environment_t* env, char* stepname, image3D* img)
{
  char filename[512];
  sprintf (filename, "./images/%s/%s.dat", env->patient_name, stepname);
  char command[1024];
  if (img->image_type == ALNSB_IMAGE_REAL)
    sprintf (command, "./scripts/display-images.sh %s %d %d %d", filename,
	     img->slices, img->rows, img->cols);
  else
    sprintf (command, "./scripts/display-images-bin.sh %s %d %d %d",
	     filename,
	     img->slices, img->rows, img->cols);
  system (command);
}

// Use '1' for the size_z for 2D images, and 1 for size_z and size_x
// for 1D images.
// Use '0' for either size_z or size_x or size_z if the size in this
// dimension is unknown.
//
// Ex: read image of 0 x 42 x 51 => read a 3D image with unknown number
// of slices, each of known size 42 x 51.
//
// Ex: read image of 1 x 0 x 51 => read a 2D image with unknown number
// of rows, each of known size 51.
//
// Ex: read image of 1 x 1 x 0 => read a 1D image with unknown number
// of elements.
//
static
image3D* load_image (s_alnsb_environment_t* env, char* stepname, int type,
		     size_t size_z, size_t size_x, size_t size_y)
{
  if (env->verbose_level > 0)
    fprintf (stderr, "[%s] Loading result from file...\n", stepname);
  char filename[512];
  sprintf (filename, "./images/%s/%s.dat", env->patient_name, stepname);
  void* data = NULL;
  size_t reads = 0;
  size_t sz = env->num_slices * env->slice_size_x * env->slice_size_y;
  size_t elt_sz = 0;
  switch (type)
    {
    case ALNSB_IMAGE_REAL: elt_sz = sizeof(ALNSB_IMAGE_TYPE_REAL); break;
    case ALNSB_IMAGE_BINARY: elt_sz = sizeof(ALNSB_IMAGE_TYPE_BIN); break;
    case ALNSB_IMAGE_INTEGER: elt_sz = sizeof(ALNSB_IMAGE_TYPE_INT); break;
    default:
      {
	fprintf (stderr, "[ERROR] Unsupported image type for %s\n", filename);
	exit (1);
      }
    }
  data = alnsb_read_data_from_binary_file_nosz (filename, elt_sz, sz, &reads);

  image3D* ret = (image3D*) malloc (sizeof(image3D));
  ret->image_type = type;
  ret->image3D = ret;
  ret->slices = size_z;
  ret->rows = size_x;
  ret->cols = size_y;
  ret->num_pixels = reads;
  ret->data = data;
  ret->object_size = reads * elt_sz;

  // Some dimension was unknown.
  if (ret->slices == 0 || ret->rows == 0 || ret->cols == 0)
    {
      // 3D image whose number of slices was unknown, e.g.
      // of size ?? x some_size x some_size.
      if (ret->slices == 0)
	ret->slices = reads / (ret->rows * ret->cols);
      // 2D image whose number of rows was unknown, e.g. 1 x ?? x some_size.
      else if (ret->slices == 1 && ret->rows == 0)
	ret->rows = reads / ret->cols;
      // 1D image whose number of elements was unknown, e.g. 1 x 1 x ??.
      else if (ret->slices == 1 && ret->rows == 1 && ret->cols == 0)
	ret->cols = reads;
      else
	{
	  fprintf (stderr,
		   "[ERROR] Unsupported size description in loader for %s\n",
		   filename);
	  exit (1);
	}
    }
  return ret;
}


static
void dump_image (s_alnsb_environment_t* env, char* stepname, image3D* img)
{
  if (env->dump_images)
    {
      char filename[512];
      sprintf (filename, "./images/%s/%s.dat", env->patient_name, stepname);
      alnsb_save_data_to_file (img->data, filename, sizeof(char),
			       img->object_size);
    }
}

static
void pass_starts (s_alnsb_environment_t* env, int pass_id)
{
  char* pass_name = env->pass_options[pass_id].pass_name;
  if (env->verbose_level)
    fprintf (stdout, "[INFO] Starting %s pass\n", pass_name);

  if (env->timer)
    alnsb_timer_start ();
}

static
void pass_ends (s_alnsb_environment_t* env, int pass_id, image3D* img)
{
  char* pass_name = env->pass_options[pass_id].pass_name;

  if (env->timer)
    {
      alnsb_timer_stop ();
      alnsb_timer_print (stdout, pass_name);
    }

  dump_image (env, pass_name, img);
  if (env->pass_options[pass_id].display_pass_result)
    display_image (env, pass_name, img);

  if (env->verbose_level)
    fprintf (stdout, "[INFO] Done with %s pass\n", pass_name);

}



// output 0: result of emtv (currently: loaded image).
void emtv_wrapper (s_alnsb_environment_t* env,
		   s_alnsb_step_t* step_data)
{
  image3DReal* output = NULL;

  int pass_id = EMTV_PASS;
  char* pass_name = env->pass_options[pass_id].pass_name;

  pass_starts (env, pass_id);

  // Step is in charge of allocating output data structure, and works
  // on concrete image types (e.g., image3DReal).
  output = (image3DReal*)
    load_image (env, pass_name, ALNSB_IMAGE_REAL,
		env->num_slices, env->slice_size_x, env->slice_size_y);

  pass_ends (env, pass_id, output->image3D);

  // Register the output in the step I/O description.
  alnsb_step_push_output (&step_data, output->image3D);
}




// input 0: result of emtv.
// output 0: result of rotation.
void rotation_wrapper (s_alnsb_environment_t* env,
		       s_alnsb_step_t* step_data)
{
  // Get the input image(s) from the step I/O description.
  image3DReal* input = (image3DReal*)step_data->read[0];
  image3DReal* output = NULL;

  int pass_id = ROTATION_PASS;
  char* pass_name = env->pass_options[pass_id].pass_name;
  int load_output = env->pass_options[pass_id].load_pass_result;

  pass_starts (env, pass_id);

  // Step is in charge of allocating output data structure, and works
  // on concrete image types (e.g., image3DReal).
  if (! load_output)
    rotation_cpu (env, input, &output);
  else
    output = (image3DReal*) load_image (env, pass_name, ALNSB_IMAGE_REAL,
					input->slices, input->rows, input->cols);

  pass_ends (env, pass_id, output->image3D);

  // Register the output in the step I/O description.
  alnsb_step_push_output (&step_data, output->image3D);
}


// input 0: result of rotation.
// output 0: result of preparation.
void levelscale_wrapper (s_alnsb_environment_t* env,
			 s_alnsb_step_t* step_data)
{
  // Get the input image(s) from the step I/O description.
  image3DReal* input = (image3DReal*)step_data->read[0];
  image3DReal* output = NULL;

  int pass_id = LEVELSCALE_PASS;
  char* pass_name = env->pass_options[pass_id].pass_name;
  int load_output = env->pass_options[pass_id].load_pass_result;

  pass_starts (env, pass_id);

  // Step is in charge of allocating output data structure, and works
  // on concrete image types (e.g., image3DReal).
  if (! load_output)
    levelscale_cpu (env, input, &output);
  else
    output = (image3DReal*) load_image (env, pass_name, ALNSB_IMAGE_REAL,
					input->slices, input->rows, input->cols);

  pass_ends (env, pass_id, output->image3D);

  // Register the output in the step I/O description.
  alnsb_step_push_output (&step_data, output->image3D);
}


// input 0: result of emtv.
// output 0: result of segmentation.
void segmentation_wrapper (s_alnsb_environment_t* env,
			   s_alnsb_step_t* step_data)
{
  // Get the input image(s) from the step I/O description.
  image3DReal* input = (image3DReal*)step_data->read[0];
  image3DReal* output = NULL;

  int pass_id = SEGMENTATION_PASS;
  char* pass_name = env->pass_options[pass_id].pass_name;
  int load_output = env->pass_options[pass_id].load_pass_result;

  pass_starts (env, pass_id);

  // Step is in charge of allocating output data structure, and works
  // on concrete image types (e.g., image3DReal).
  if (! load_output)
    segmentation_cpu (env, input, &output);
  else
    output = (image3DReal*) load_image (env, pass_name, ALNSB_IMAGE_REAL,
					input->slices, input->rows, input->cols);

  pass_ends (env, pass_id, output->image3D);

  // Register the output in the step I/O description.
  alnsb_step_push_output (&step_data, output->image3D);
}


// input 0: result of segmentation.
// output 0: result of segmentation mask.
void segmentationMask_wrapper (s_alnsb_environment_t* env,
			       s_alnsb_step_t* step_data)
{
  // Get the input image(s) from the step I/O description.
  image3DReal* input = (image3DReal*)step_data->read[0];
  image3DBin* output = NULL;

  int pass_id = SEGMENTATIONMASK_PASS;
  char* pass_name = env->pass_options[pass_id].pass_name;
  int load_output = env->pass_options[pass_id].load_pass_result;

  pass_starts (env, pass_id);

  // Step is in charge of allocating output data structure, and works
  // on concrete image types (e.g., image3DReal).
  if (! load_output)
    segmentationMask_cpu (env, input, &output);
  else
    output = (image3DBin*) load_image (env, pass_name, ALNSB_IMAGE_BINARY,
				       input->slices, input->rows, input->cols);

  pass_ends (env, pass_id, output->image3D);

  // Register the output in the step I/O description.
  alnsb_step_push_output (&step_data, output->image3D);
}



// input 0: result of segmentationMask.
// output 0: result of preselection.
void preselection_wrapper (s_alnsb_environment_t* env,
			   s_alnsb_step_t* step_data)
{
  // Get the input image(s) from the step I/O description.
  image3DBin* input = (image3DBin*)step_data->read[0];
  image3DBin* output = NULL;

  int pass_id = PRESELECTION_PASS;
  char* pass_name = env->pass_options[pass_id].pass_name;
  int load_output = env->pass_options[pass_id].load_pass_result;

  pass_starts (env, pass_id);

  // Step is in charge of allocating output data structure, and works
  // on concrete image types (e.g., image3DReal).
  if (! load_output)
    preselection_cpu (env, input, &output);
  else
    output = (image3DBin*) load_image (env, pass_name, ALNSB_IMAGE_BINARY,
				       input->slices, input->rows, input->cols);

  pass_ends (env, pass_id, output->image3D);

  // Register the output in the step I/O description.
  alnsb_step_push_output (&step_data, output->image3D);
}


// input 0: result of preselection.
// output 0: result of featureExtraction.
void featureExtraction_wrapper (s_alnsb_environment_t* env,
				s_alnsb_step_t* step_data)
{
  // Get the input image(s) from the step I/O description.
  image3DReal* inputPrep = (image3DReal*)step_data->read[0];
  image3DBin* inputPres = (image3DBin*)step_data->read[1];
  image3DReal* output = NULL;
  fprintf (stdout, "output pointer initilize\n");

  int pass_id = FEATUREEXTRACTION_PASS;
  char* pass_name = env->pass_options[pass_id].pass_name;
  int load_output = env->pass_options[pass_id].load_pass_result;

  pass_starts (env, pass_id);

  // Step is in charge of allocating output data structure, and works
  // on concrete image types (e.g., image3DReal).
  if (! load_output)
    featureExtraction_cpu (env, inputPrep, inputPres, &output);
  else
    output = (image3DReal*) load_image (env, pass_name, ALNSB_IMAGE_REAL,
					1, 0, env->classifier_num_features);

  pass_ends (env, pass_id, output->image3D);

  // Register the output in the step I/O description.
  alnsb_step_push_output (&step_data, output->image3D);
}



// input 0: result of preselection.
// input 1: result of featureExtraction.
// output 0: result of preselection.
void classification_wrapper (s_alnsb_environment_t* env,
			     s_alnsb_step_t* step_data)
{
  // Get the input image(s) from the step I/O description.
  image3DReal* inputPrep = (image3DReal*)step_data->read[0];
  image3DBin* inputPres = (image3DBin*)step_data->read[1];
  image3DReal* features = (image3DReal*)step_data->read[2];
  image3DReal* output = NULL;

  int pass_id = CLASSIFICATION_PASS;
  char* pass_name = env->pass_options[pass_id].pass_name;
  int load_output = env->pass_options[pass_id].load_pass_result;

  pass_starts (env, pass_id);

  // Step is in charge of allocating output data structure, and works
  // on concrete image types (e.g., image3DReal).
  if (! load_output)
    classification_cpu (env, inputPrep, inputPres, features, &output);
  else
    output = (image3DReal*) load_image (env, pass_name, ALNSB_IMAGE_REAL,
					inputPrep->slices,
					inputPrep->rows, inputPrep->cols);

  pass_ends (env, pass_id, output->image3D);

  // Register the output in the step I/O description.
  alnsb_step_push_output (&step_data, output->image3D);
}




void alnsb_pipeline (s_alnsb_environment_t* env)
{
  s_alnsb_step_t* emtv_io = NULL;
  s_alnsb_step_t* rot_io = NULL;
  s_alnsb_step_t* prep_io = NULL;
  s_alnsb_step_t* seg_io = NULL;
  s_alnsb_step_t* segmask_io = NULL;
  s_alnsb_step_t* presel_io = NULL;
  s_alnsb_step_t* featExt_io = NULL;
  s_alnsb_step_t* class_io = NULL;
  size_t i;

  fprintf (stdout, "* * * * * * * * * CAD Pipeline starts * * * * * * * * * *\n");
  
  // env::(emtv)
  // (fake input to "start" the graph)
  alnsb_step_push_input (&emtv_io, NULL);

  // emtv -> [emtv_io]
  emtv_wrapper (env, emtv_io);
  alnsb_step_push_input (&rot_io, emtv_io->write[0]);

  // Stage 0: preparation.
  // [emtv_io] -> (rotation) -> [rot_io]
  rotation_wrapper (env, rot_io);
  alnsb_step_push_input (&prep_io, rot_io->write[0]);

  // [rot_io] -> (levelscale) -> [prep_io]
  levelscale_wrapper (env, prep_io);
  alnsb_step_push_input (&seg_io, prep_io->write[0]);
  alnsb_step_push_input (&featExt_io, prep_io->write[0]);
  alnsb_step_push_input (&class_io, prep_io->write[0]);

  // Stage 1: segmentation.
  // [prep_io] -> (segmentation) -> [segmask_io];
  segmentation_wrapper (env, seg_io);
  alnsb_step_push_input (&segmask_io, seg_io->write[0]);


  // Stage 2: build segmentation mask.
  // [segmask_io] -> (segmentationMask) -> [segmasked_io];
  segmentationMask_wrapper (env, segmask_io);
  alnsb_step_push_input (&presel_io, segmask_io->write[0]);


  // Stage 3: preselection.
  // [segmasked_io] -> (preselection) -> [presel_io];
  preselection_wrapper (env, presel_io);
  alnsb_step_push_input (&class_io, presel_io->write[0]);
  alnsb_step_push_input (&featExt_io, presel_io->write[0]);


  // Stage 4: feature extraction.
  // [rot_io],[presel_io] -> (featureExtraction) -> [features_io];
  featureExtraction_wrapper (env, featExt_io);
  alnsb_step_push_input (&class_io, featExt_io->write[0]);

  //fprintf (stdout, "--segmentation             ==>      92.625 seconds\n");
  //fprintf (stdout, "--segmentationMask         ==>      184.564 seconds\n");
  //fprintf (stdout, "--preslection              ==>      20.122 seconds\n");
  //fprintf (stdout, "--featureExtraction        ==>      18.867 seconds\n");
  //fprintf (stdout, "--classification           ==>      5.344 seconds\n");

  // Stage 5: classification.
  // [rot_io],[presel_io],[features_io] -> (classification) -> [output_io];
  classification_wrapper (env, class_io);


  fprintf (stdout, "* * * * * * * * * CAD Pipeline ends * * * * * * * * * *\n");

  // Be clean.
  alnsb_step_free (emtv_io);
  alnsb_step_free (rot_io);
  alnsb_step_free (prep_io);
  alnsb_step_free (seg_io);
  alnsb_step_free (segmask_io);
  alnsb_step_free (presel_io);
  alnsb_step_free (featExt_io);
  alnsb_step_free (class_io);
}
