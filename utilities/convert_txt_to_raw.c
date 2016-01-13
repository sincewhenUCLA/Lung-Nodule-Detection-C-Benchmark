/**
 * convert_txt_to_raw.c: this file is part of the ALNSB project.
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
void save_data_to_file (void* data, char* filename, unsigned int elt_sz,
			unsigned int elt_count)
{
  FILE* f = fopen (filename, "w+");
  if (f == NULL)
    {
      fprintf (stderr, "[ERROR] impossible to create file (system full?)\n");
      exit (1);
    }
  fwrite (data, elt_sz, elt_count, f);
  fclose (f);
}




/**
 * Convert an ASCII float image to raw data.
 *
 */
int main(int argc, char** argv)
{
  if (argc != 4)
    {
      printf ("Usage: %s <input_file.txt> <output_file.dat> <nb_pix>\n",
	      argv[0]);
      printf ("=> Converts a txt file containing ASCII floats into raw data\n");
      exit (1);
    }
  unsigned int nb_pix = atoi (argv[3]);
  printf ("[ASCIIToRaw] Converting %s (txt) of %d elements to %s (raw/float)...\n",
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
  float* data_p = data;
  char buffer[32];
  char* tmp;
  unsigned int i;
  int nc = 0;
  char b[2];
  for (i = 0; i < nb_pix; ++i)
    {
      int pos = 0;
      while ((nc = fread (b, sizeof(char), 1, fr)))
	{
	  if ((b[0] >= '0' && b[0] <= '9') || b[0] == '.')
	    buffer[pos++] = b[0];
	  else
	    break;
	}
      if (nc == 0)
	break;
      if (pos == 0)
	{
	  --i;
	  continue;
	}
      buffer[pos] = '\0';
      float val = atof (buffer);
      *(data_p++) = val;
    }


  int remainder = 0;
  while ((nc = fread (buffer, sizeof(char), 32, fr)))
    remainder += nc;

  fclose (fr);

  if (remainder > 0)
    printf ("[WARNING] remainder: %d characters not read from the file!!!\n",
	    remainder);

  if (i != nb_pix)
    {
      printf ("[ERROR] read %d pixels, expected %d\n", i, nb_pix);
      exit (1);
    }

  save_data_to_file (data, argv[2], sizeof(float), nb_pix);

  free (data);

  return 0;
}
