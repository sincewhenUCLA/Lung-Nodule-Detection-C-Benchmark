## filename: Makefile.am
##
##
## Copyright (c) 2014, 2015 University of California Los Angeles. All
## rights reserved.
##
## Written by Shiwen Shen, Prashant Rawat and Louis-Noel Pouchet
##

## Define the compiler binaries
CC=gcc-mp-4.9
CXX=g++-mp-4.9
## Define the compilation flags for optimized and debug cases.
INCLUDES=-I. -Iutilities -Itoolbox
CFLAGS_OPT=-O3 -std=c99 -lm  $(INCLUDES)  -D__ALNSB_RESTRICT_PTR=restrict
CFLAGS_OPT_OMP=-O3 -std=c99 -lm -fopenmp $(INCLUDES)  -D__ALNSB_RESTRICT_PTR=restrict
CFLAGS_DEBUG= -std=c99 -lm -g -ggdb $(INCLUDES) -D__ALNSB_RESTRICT_PTR=restrict
## Select which compilation mode (optimized, parallel or debug)
#CFLAGS=$(CFLAGS_OPT)
CFLAGS=$(CFLAGS_OPT_OMP)
## Output binary name for the reference pipeline.
PROG_NAME = alnsb

UTILITIES_SRC =					\
	utilities/environment.c			\
	utilities/images.c			\
	utilities/memfuncs.c			\
	utilities/step.c			\
	utilities/timer.c			\
	utilities/file_io.c

TOOLBOX_SRC =					\
	toolbox/rotate.c			\
	toolbox/level.c				\
	toolbox/scale.c				\
	toolbox/erode.c				\
	toolbox/dilate.c			\
	toolbox/floodfill.c			\
	toolbox/bwconncomp.c			\
	toolbox/imPerimeter.c			\
	toolbox/imSurface.c			\
	toolbox/imMeanBreadth.c			\
	toolbox/imEuler3d.c			\
	toolbox/stdev.c				\
	toolbox/skewness.c			\
	toolbox/kurtosis.c

STAGES_SRC =							\
	stages/rotation/rotation_step.c				\
	stages/levelscale/levelscale_step.c			\
	stages/segmentation/segmentation_step.c			\
	stages/segmentationMask/segmentationMask_step.c		\
	stages/preselection/preselection_step.c			\
	stages/featureExtraction/featureExtraction_step.c	\
	stages/classification/classification_step.c


DRIVER_SRC =					\
	driver/main.c				\
	driver/pipeline.c			\
	driver/options.c


OBJECTS_BASE = $(UTILITIES_SRC:.c=.o) $(DRIVER_SRC:.c=.o) $(TOOLBOX_SRC:.c=.o) $(STAGES_SRC:.c=.o)


all: $(PROG_NAME) convert_txt_to_raw convert_raw_to_txt

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@
.cc.o:
	$(CXX) $(CFLAGS) -c $< -o $@

$(PROG_NAME): $(OBJECTS_BASE)
	$(CC) $(CFLAGS) $(OBJECTS_BASE) -o $(PROG_NAME)

clean:
	rm -f $(PROG_NAME) $(OBJECTS_BASE) convert_txt_to_raw convert_raw_to_txt

convert_txt_to_raw: utilities/convert_txt_to_raw.c
	$(CC) $(CFLAGS) utilities/convert_txt_to_raw.c -o convert_txt_to_raw

convert_raw_to_txt: utilities/convert_raw_to_txt.c
	$(CC) $(CFLAGS) utilities/convert_raw_to_txt.c -o convert_raw_to_txt

unzip-images: 
	cd images/NLST_R0960B_OUT4 && tar xzf emtv.tar.gz

run: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 1 --thresolding_skip4 --rotate-90-right --enable-timing --classification-display


run-feat: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 42 --thresolding_skip4 --rotate-90-right --enable-timing --classification-display --rotation-skip --levelscale-skip --segmentation-skip --segmentationMask-skip --preselection-skip

run-class: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 42 --thresolding_skip4 --rotate-90-right --enable-timing --classification-display --rotation-skip --levelscale-skip --segmentation-skip --segmentationMask-skip --preselection-skip --featureExtraction-skip

run-test: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 42 --thresolding_skip4 --rotate-90-right --enable-timing --rotation-skip --levelscale-skip --segmentation-skip --segmentationMask-skip --preselection-skip --featureExtraction-skip --classification-skip

run-segmask: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 42 --thresolding_skip4 --rotate-90-right --enable-timing --rotation-skip --levelscale-skip --segmentation-skip

run-test-convert: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 42 --patient_folder NLST_test_converter --num_slices 124 --slice-size-x 716 --slice-size-y 716 --enable-timing --emtv-display --levelscale-display --segmentationMask-display --preselection-display

run-test-lidc: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 42 --patient_folder LIDCtest --num_slices 133 --slice-size-x 512 --slice-size-y 512 --enable-timing --emtv-display --levelscale-display --segmentationMask-display --preselection-display --segmentation_ulab0 0.3 --segmentation_ulab1 0.6

run-lidc-16: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 1 --thresold_upperBand 1024 --thresold_lowerBand -1024 --patient_folder LIDC-IDRI-0016 --num_slices 141 --slice-size-x 512 --slice-size-y 512 --enable-timing --levelscale-display --segmentationMask-display --segmentation_ulab0 0.3 --segmentation_ulab1 0.6 > lidc_16_log.txt

run-lidc-31: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 1 --patient_folder LIDC-IDRI-0031 --num_slices 133 --slice-size-x 512 --slice-size-y 512 --enable-timing --levelscale-display --segmentationMask-display --segmentation_ulab0 0.3 --segmentation_ulab1 0.6 > lidc_31_log.txt

run-lidc-113: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 1 --patient_folder LIDC-IDRI-0113 --num_slices 145 --slice-size-x 512 --slice-size-y 512 --enable-timing --levelscale-display --segmentationMask-display --segmentation_ulab0 0.3 --segmentation_ulab1 0.6 > lidc_113_log.txt

run-lidc-116: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 1 --patient_folder LIDC-IDRI-0116 --num_slices 133 --slice-size-x 512 --slice-size-y 512 --enable-timing --levelscale-display --segmentationMask-display --segmentation_ulab0 0.3 --segmentation_ulab1 0.6 > lidc_116_log.txt

run-lidc-122: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 1 --patient_folder LIDC-IDRI-0122 --num_slices 119 --slice-size-x 512 --slice-size-y 512 --enable-timing --levelscale-display --segmentationMask-display --segmentation_ulab0 0.3 --segmentation_ulab1 0.6 > lidc_122_log.txt

run-lidc-126: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 1 --patient_folder LIDC-IDRI-0126 --num_slices 137 --slice-size-x 512 --slice-size-y 512 --enable-timing --levelscale-display --segmentationMask-display --segmentation_ulab0 0.3 --segmentation_ulab1 0.6 > lidc_126_log.txt

run-lidc-145: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 1 --patient_folder LIDC-IDRI-0145 --num_slices 139 --slice-size-x 512 --slice-size-y 512 --enable-timing --levelscale-display --segmentationMask-display --segmentation_ulab0 0.3 --segmentation_ulab1 0.6 > lidc_145_log.txt

run-lidc-146: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 1 --patient_folder LIDC-IDRI-0146 --num_slices 310 --slice-size-x 512 --slice-size-y 512 --enable-timing --levelscale-display --segmentationMask-display --segmentation_ulab0 0.3 --segmentation_ulab1 0.6 > lidc_146_log.txt

run-lidc-155: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 1 --patient_folder LIDC-IDRI-0155 --num_slices 238 --slice-size-x 512 --slice-size-y 512 --enable-timing --levelscale-display --segmentationMask-display --segmentation_ulab0 0.3 --segmentation_ulab1 0.6 > lidc_155_log.txt

run-lidc-157: $(PROG_NAME)
	./alnsb --show-environment --verbose-level 1 --patient_folder LIDC-IDRI-0157 --num_slices 133 --slice-size-x 512 --slice-size-y 512 --enable-timing --levelscale-display --segmentationMask-display --segmentation_ulab0 0.3 --segmentation_ulab1 0.6 > lidc_157_log.txt












