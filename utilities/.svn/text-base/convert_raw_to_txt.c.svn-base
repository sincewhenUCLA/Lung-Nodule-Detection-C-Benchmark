/**
 * convert_raw_to_txt.c: this file is part of the ALNSB project.
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
#include <string.h>

static
void save_data_to_file_ascii (float* data, char* filename,
			      unsigned int elt_count)
{
  FILE* f = fopen (filename, "w+");
  if (f == NULL)
    {
      fprintf (stderr, "[ERROR] impossible to create file (system full?)\n");
      exit (1);
    }
  unsigned int i;
  char buffer[64];
  for (i = 0; i < elt_count; ++i)
    {
      sprintf (buffer, "%f ", *(data++));
      fwrite (buffer, sizeof(char), strlen(buffer), f);
    }
  fclose (f);
}




/**
 * Convert float image to an ASCII image.
 *
 */
int main(int argc, char** argv)
{
  if (argc != 4)
    {
      printf ("Usage: %s <input_file.dat> <output_file.txt> <nb_pix>\n",
	      argv[0]);
      printf ("=> Converts a txt file containing raw data into ASCII floats\n");
      exit (1);
    }

  unsigned int nb_pix = atoi (argv[3]);
  printf ("[RawToASCII] Converting %s (raw/float) of %d elements to %s (txt)...\n",
	  argv[1], nb_pix, argv[2]);

  FILE* fr = fopen(argv[1], "r");
  if (fr == NULL)
    {
      printf ("[ERROR] File %s cannot be opened\n", argv[1]);
      exit (1);
    }
  if (nb_pix == 0)
    {
      printf ("[ERROR] Nb pixels cannot be 0\n");
      exit (1);
    }

  float* data = (float*) malloc (sizeof(float) * nb_pix);
  int nb_elts_read = fread (data, sizeof(float), nb_pix, fr);
  fclose (fr);

  if (nb_elts_read != nb_pix)
    {
      printf ("[ERROR] read %d pixels, expected %d\n", nb_elts_read, nb_pix);
      exit (1);
    }

  save_data_to_file_ascii (data, argv[2], nb_pix);

  free (data);

  return 0;
}
