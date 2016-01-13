/**
 * file_io.c: this file is part of the ALNSB project.
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
#include <utilities/file_io.h>

static
void* alnsb_read_data_from_binary_file_ (char* filename, size_t elt_sz,
					 size_t elt_count, int noerror,
					 size_t* read_count)
{
  FILE* f = fopen (filename, "r");
  if (f == NULL)
    {
      fprintf (stderr, "[ERROR] impossible to open file %s\n", filename);
      exit (1);
    }
  void* ret = malloc (elt_sz * elt_count);
  if (ret == NULL)
    {
      fprintf (stderr, "[ERROR] Memory exhausted\n");
      exit (1);
    }
  char* cur = ret;
#define BUFF_SZ 1024
  char buffer[BUFF_SZ];
  int nc;
  size_t nb_read = 0;
  while ((nc = fread (buffer, sizeof(char), BUFF_SZ, f)))
    {
      int i;
      for (i = 0; i < nc; ++i)
	{
	  if (nb_read < elt_sz * elt_count)
	    *(cur++) = buffer[i];
	  ++nb_read;
	}
    }
  fclose (f);

  if (elt_sz * elt_count != nb_read && !noerror)
    fprintf (stderr, "[WARNING] Loaded %d characters from a file containing %d characters\n", elt_sz * elt_count, nb_read);
  if (elt_sz * elt_count != nb_read)
    ret = realloc (ret, nb_read);
  if (read_count)
    *read_count = nb_read / elt_sz;

  return ret;
}

// Standard reader when the number of elements is known.
void* alnsb_read_data_from_binary_file (char* filename, size_t elt_sz,
					size_t elt_count)
{
  return alnsb_read_data_from_binary_file_ (filename, elt_sz,
					    elt_count, 0, NULL);
}

// Use when number of elements is unknown, giving a max. elt_count possible.
void* alnsb_read_data_from_binary_file_nosz (char* filename,
					     size_t elt_sz,
					     size_t max_elt_count,
					     size_t* size)
{
  return alnsb_read_data_from_binary_file_ (filename, elt_sz,
					    max_elt_count, 1, size);
}

// Standard saver to binary format.
void alnsb_save_data_to_file (void* data, char* filename, size_t elt_sz,
			      size_t elt_count)
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

// Standard saver to text format.
void alnsb_save_data_to_file_ascii (float* data, char* filename,
				    size_t elt_count)
{
  FILE* f = fopen (filename, "w+");
  if (f == NULL)
    {
      fprintf (stderr, "[ERROR] impossible to create file (system full?)\n");
      exit (1);
    }
  size_t i;
  char buffer[32];
  for (i = 0; i < elt_count; ++i)
    {
      sprintf (buffer, "%f ", *(data++));
      fwrite (buffer, sizeof(float), strlen(buffer), f);
    }
  fclose (f);
}
