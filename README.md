# Lung-Nodule-Detection-C-Benchmark

=======================================================
* ALNSB: the Adaptive Lung Nodule Screening Benchmark *
=======================================================

v 0.1.

Contact: Alex Bui <buia@mii.ucla.edu>

Written by: Shiwen Shen, Prashant Rawat, Louis-Noel Pouchet and William Hsu.



* Instructions:
---------------

1) Edit Makefile and update the CC and CXX variables to point to a
working GCC version. Also edit scripts/display-image.sh and
scripts/display-images-bin.sh to use a valid command line to display
an image from raw floating point data (32 bits). Default setup is for
Mac OS X using 'convert' (imageMagick) and 'open' to open a pdf file.

2) Run, for the first-time install or if the emtv input image has been
damaged:

$> make unzip-images

This will expand the reference emtv output for patient
NLST_R0960B_OUT4 in the images/NLST_R0960B_OUT4 directory. All images
must be stored in the directory 'images/<patient>' where <patient> is
a unique name (e.g., NLST_R0960B_OUT4) and can be specified using the
--patient-folder option of the program.

3) To run the full pipeline:

$> make run

4) To change algorithm parameters, behavior, etc. run the binary directly:

To get available options:
$> ./alnsb --help

Then, for example:
$> ./alnsb --show-environment --classification-display

5) To convert a text image/matrix (e.g., a sequence of floating point
numbers separated by spaces/newlines in plain ascii format) to a valid
input to the pipeline (raw/float format), do:

$> ./convert_txt_to_raw input.txt output.dat 42

where 42 is the number of distinct numbers in the text file.

For the reverse, that is converting an input/output of the pipeline
into plain ascii text format, do:
$> ./convert_raw_to_ascii input.dat output.txt 42


These programs can be used seamlessly to convert binary or real
images, as well as the matrices input to the classifier (they must be
stored in float/raw format too, like images, however can be stored
anywhere).

6) To display an image of 124 slices of size 716x716, in float/raw
format:

$> scripts/display-image.sh image.dat 124 716 716

(beware to edit this script first as needed, as said above).


