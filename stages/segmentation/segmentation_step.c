/**
 * segmentation_step.c: this file is part of the ALNSB project.
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
#include <math.h>
#include <stages/segmentation/segmentation_step.h>

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))


void segmentation_cpu (s_alnsb_environment_t* __ALNSB_RESTRICT_PTR env,
		       image3DReal* __ALNSB_RESTRICT_PTR input,
		       image3DReal** __ALNSB_RESTRICT_PTR output)
{
  // Allocate output data.
  *output = image3DReal_alloc (input->slices, input->rows, input->cols);

  // Set up algorithm parameters set by the user.
  float lp = env->segmentation_lp;
  float errb[] = {env->segmentation_errb[0], env->segmentation_errb[1]};
  float ulab[] = {env->segmentation_ulab[0], env->segmentation_ulab[1]};
  float cc = env->segmentation_cc;
  float c_convergence = env->segmentation_c_convergence;
  float steps = env->segmentation_steps;
  float beta = env->segmentation_beta;
  float err_c = env->segmentation_err_c;
  float max_err = env->segmentation_max_err;


  // Locals.
  int slices = input->slices;
  int rows = input->rows;
  int cols = input->cols;
  int t, i, j, k, iter = 0;
  float ulab_p[2];
  float alpha = lp;
  // Set to 1 for printing some debug info while executing.
  int debug = 0;

  // 3D view of input/output arrays.
  ALNSB_IMRealTo3D(input,ur);
  ALNSB_IMRealTo3D(*output,ut);

  // Allocate local 3D arrays.
  ALNSB_REAL3D_alloc(Cs, slices, rows, cols);
  ALNSB_REAL3D_alloc(Ct, slices, rows, cols);
  ALNSB_REAL3D_alloc(pp1, slices, rows, cols+1);
  ALNSB_REAL3D_alloc(pp2, slices, rows+1, cols);
  ALNSB_REAL3D_alloc(pp3, slices+1, rows, cols);
  ALNSB_REAL3D_alloc(u, slices, rows, cols);
  ALNSB_REAL3D_alloc(divp, slices, rows, cols);
  ALNSB_REAL3D_alloc(ps, slices, rows, cols);
  ALNSB_REAL3D_alloc(pt, slices, rows, cols);
  ALNSB_REAL3D_alloc(gk, slices, rows, cols);
  ALNSB_REAL3D_alloc(pts, slices, rows, cols);
  ALNSB_REAL3D_alloc(erru, slices, rows, cols);

#pragma omp parallel for private(j,k)
  for (i = 0; i < slices; i++) {
     for (j = 0; j < rows; j++) {
        for (k = 0; k < cols; k++) {
          Cs[i][j][k] = pow (fabs (ur[i][j][k] - ulab[0]), beta);
          Ct[i][j][k] = pow (fabs (ur[i][j][k] - ulab[1]), beta);
          u[i][j][k] = (Cs[i][j][k] - Ct[i][j][k]) >= 0;
          ps[i][j][k] = pt[i][j][k] = min (Cs[i][j][k], Ct[i][j][k]);
        }
     }
  }

  while (err_c > errb[1]) {
      if (err_c < c_convergence)
         errb[0] = errb[1];

      iter++;
      if (debug)
	printf ("outer loop iteration no: %d\n", iter);

      for (t = 0; t < env->segmentation_max_steps; t++) {
	if (debug)
	    printf ("inner loop iteration no: %d\n", t+1);
#pragma omp parallel for private(j,k)
          for (i = 0; i < slices; i++)
            for (j = 0; j < rows; j++)
              for (k = 0; k < cols; k++)
                pts[i][j][k] = divp[i][j][k] - (ps[i][j][k] - pt[i][j][k] + (u[i][j][k] / cc));

#pragma omp parallel for private(j,k)
          for (i = 0; i < slices; i++)
            for (j = 0; j < rows; j++)
              for (k = 1; k < cols; k++)
                pp1[i][j][k] = pp1[i][j][k] + (steps * (pts[i][j][k] - pts[i][j][k-1]));

#pragma omp parallel for private(j,k)
          for (i = 0; i < slices; i++)
            for (j = 1; j < rows; j++)
              for (k = 0; k < cols; k++)
                pp2[i][j][k] = pp2[i][j][k] + (steps * (pts[i][j][k] - pts[i][j-1][k]));

#pragma omp parallel for private(j,k)
          for (i = 1; i < slices; i++)
            for (j = 0; j < rows; j++)
              for (k = 0; k < cols; k++)
                pp3[i][j][k] = pp3[i][j][k] + (steps * (pts[i][j][k] - pts[i-1][j][k]));

#define alnsb_power2(a) ((a)*(a))

#pragma omp parallel for private(j,k)
          for (i = 0; i < slices-1; i++) {
            for (j = 0; j < rows-1; j++) {
              for (k = 0; k < cols-1; k++) {
		// power of 2 was hardwired, beta was not used here.
                gk[i][j][k] = sqrt ((alnsb_power2(pp1[i][j][k]) + alnsb_power2(pp1[i][j][k+1]) + alnsb_power2(pp2[i][j][k]) + alnsb_power2(pp2[i][j+1][k]) + alnsb_power2(pp3[i][j][k]) + alnsb_power2(pp3[i+1][j][k])) * 0.5);
                gk[i][j][k] = (gk[i][j][k] <= alpha) + (gk[i][j][k] > alpha) * (gk[i][j][k] / alpha);
                gk[i][j][k] = 1 / gk[i][j][k];
              }
            }
          }

#pragma omp parallel for private(j,k)
          for (i = 0; i < slices; i++)
            for (j = 0; j < rows; j++)
              for (k = 1; k < cols; k++)
                pp1[i][j][k] = (0.5 * (gk[i][j][k] + gk[i][j][k-1])) * pp1[i][j][k];

#pragma omp parallel for private(j,k)
          for (i = 0; i < slices; i++)
            for (j = 1; j < rows; j++)
              for (k = 0; k < cols; k++)
                pp2[i][j][k] = (0.5 * (gk[i][j][k] + gk[i][j-1][k])) * pp2[i][j][k];

#pragma omp parallel for private(j,k)
          for (i = 1; i < slices; i++)
            for (j = 0; j < rows; j++)
              for (k = 0; k < cols; k++)
                pp3[i][j][k] = (0.5 * (gk[i][j][k] + gk[i-1][j][k])) * pp3[i][j][k];

	  max_err = 0.0;
#pragma omp parallel for private(j,k)
          for (i = 0; i < slices-1; i++) {
            for (j = 0; j < rows-1; j++) {
              for (k = 0; k < cols-1; k++) {
                divp[i][j][k] = pp2[i][j+1][k] - pp2[i][j][k] + pp1[i][j][k+1] - pp1[i][j][k] + pp3[i+1][j][k] - pp3[i][j][k];
                pts[i][j][k] = divp[i][j][k] - (u[i][j][k] / cc) + pt[i][j][k] + (1/cc);
                Cs[i][j][k] = pow (fabs (ur[i][j][k] - ulab[0]), beta);
                ps[i][j][k] = min (pts[i][j][k], Cs[i][j][k]);
                pts[i][j][k] = ps[i][j][k] + (u[i][j][k] / cc) - divp[i][j][k];
                Ct[i][j][k] = pow (fabs (ur[i][j][k] - ulab[1]), beta);
                pt[i][j][k] = min (pts[i][j][k], Ct[i][j][k]);
                erru[i][j][k] = cc * (divp[i][j][k] + pt[i][j][k] - ps[i][j][k]);
                u[i][j][k] = u[i][j][k] - erru[i][j][k];
              }
            }
          }

#pragma omp parallel for private(j,k) reduction(+:max_err)
          for (i = 0; i < slices; i++)
            for (j = 0; j < rows; j++)
              for (k = 0; k < cols; k++)
                 max_err += fabs (erru[i][j][k]);
	  if (debug)
	    printf ("max_err = %.16f\n", max_err);
          if ((max_err / (rows*cols*slices)) < errb[0]) break;
      }
#pragma omp parallel for private(j,k)
      for (i = 0; i < slices; i++)
         for (j = 0; j < rows; j++)
            for (k = 0; k < cols; k++)
		ut[i][j][k] = u[i][j][k] > 0.5;

      ulab_p[0] = ulab[0];
      ulab_p[1] = ulab[1];
      if (debug)
	printf ("ulab_p[0]=%f, ulab_p[1]=%f\n", ulab_p[0], ulab_p[1]);
      ulab[0] = ulab[1] = 0.0;
      for (i = 0; i < slices; i++) {
	 float f = 0.0, g = 0.0, m = 0.0, n = 0.0;
#pragma omp parallel for private(k) reduction(+:f,m,g,n)
         for (j = 0; j < rows; j++) {
	   for (k = 0; k < cols; k++) {
               f += ur[i][j][k] * (1 - ut[i][j][k]);
	       m += ur[i][j][k] * ut[i][j][k];
	       g += 1 - ut[i][j][k];
               n += ut[i][j][k];
            }
         }
	 if (g > 0)
	   ulab[0] += (f/g);
	 if (n > 0)
	   ulab[1] +=  (m/n);
      }
      ulab[0] = ulab[0] / slices;
      ulab[1] = ulab[1] / slices;
      err_c = fabs (ulab_p[0] - ulab[0]) + fabs (ulab_p[1] - ulab[1]);
      if (debug)
	printf ("err_c = %f\n", err_c);
  }


  // Free temporaries.
  ALNSB_REAL3D_free(Cs);
  ALNSB_REAL3D_free(Ct);
  ALNSB_REAL3D_free(pp1);
  ALNSB_REAL3D_free(pp2);
  ALNSB_REAL3D_free(pp3);
  ALNSB_REAL3D_free(u);
  ALNSB_REAL3D_free(divp);
  ALNSB_REAL3D_free(ps);
  ALNSB_REAL3D_free(pt);
  ALNSB_REAL3D_free(gk);
  ALNSB_REAL3D_free(pts);
  ALNSB_REAL3D_free(erru);
}
