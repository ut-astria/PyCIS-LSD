/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/interface.c: pipe from pycis.c to lib/lsd*.c

Benjamin Feuge-Miller: benjamin.g.miller@utexas.edu
The University of Texas at Austin, 
Oden Institute Computational Astronautical Sciences and Technologies (CAST) group
*Date of Modification: September 03, 2021

#--------------------------------------------------------------------------------------
#PyCIS-LSD: An a-contrario detection sub-algorithm for extracting narrow lines within dense optical data cubes.
#Copyright (C) 2022, Benjamin G. Feuge-Miller, <benjamin.g.miller@utexas.edu>
#
#PyCIS-LSD is free software: you can redistribute it and/or modify 
#it under the terms of the GNU General Public License as published 
#by the Free Software Foundation, either version 3 of the License, 
#or (at your option) any later version.
#
#PyCIS-LSD is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU Affero General Public License for more details.
# 
#You should have received a copy of the GNU Affero General Public License
#along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#**NOTICE**: 
#PyCIS-LSD is modified from the source code of LSDSAR:
#"LSDSAR, a Markovian a contrario framework for line segment detection in SAR images"
#by Chenguang Liu, Rémy Abergel, Yann Gousseau and Florence Tupin. 
#Pattern Recognition, 2019).
#https://doi.org/10.1016/j.patcog.2019.107034
#*Date of Modification: April 30, 2021*
#
#**NOTICE**: 
#LSDSAR is modified from the source code of LSD:
#"LSD: a Line Segment Detector" by Rafael Grompone von Gioi,
#Jeremie Jakubowicz, Jean-Michel Morel, and Gregory Randall,
#Image Processing On Line, 2012. DOI:10.5201/ipol.2012.gjmr-lsd
#http://dx.doi.org/10.5201/ipol.2012.gjmr-lsd
#*Date of Modification: 27/06/2018*
#--------------------------------------------------------------------------------------
#------------------------------------------------------------------------------*/   

/*----------------------------------------------------------------------------*/
/*---------------------------------- Import ---------------------------------*/
/*----------------------------------------------------------------------------*/



#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include <tgmath.h>
#include <limits.h>
#include <float.h>
#include<string.h>
#include <time.h>
#include <sys/mman.h>

#include "interface.h"
#include "constants.h"
#include "misc.h"
#include "tuples.h"
#include "markov.h"
#include "lsd2.h"
#include "lsd3.h"
#include "lsd3b.h"

/*----------------------------------------------------------------------------*/
/*----------------------- LSD-WRAPPER INTERFACES -----------------------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/** LSD Simple Interface.
 */
double * lsd(int * n_out, double * img, int X, int Y,double *inputv, double inputv_size)
{
  double ang_th;   /* Gradient angle tolerance in degrees.           */
  ang_th=inputv[4];
  double log_eps = 0.0;     /* Detection threshold: -log10(NFA) > log_eps     */
  log_eps=-log10(inputv[1]);
  
  double density_th = 0.0;  /* Minimal density of region points in rectangle. */
  density_th=inputv[2];
  int n_bins = 1024;        /* Number of bins in pseudo-ordering of gradient
                               modulus.                                       */

  return LineSegmentDetection( n_out, img, X, Y,
                              ang_th, log_eps, density_th, n_bins,
                              NULL,NULL,NULL ,inputv,inputv_size);
}

/*----------------------------------------------------------------------------*/
/** LSD3 Simple Interface.
 */
double * lsd3(int * n_out, double * img, int X, int Y, int Z, 
              double *inputv, double inputv_size,double * inputvorth)
{ 
  double ang_th;   /* Gradient angle tolerance in degrees.           */
  ang_th=inputv[4];
  double log_eps = 0.0;     /* Detection threshold: -log10(NFA) > log_eps     */
  log_eps=-log10(inputv[1]);

  double density_th = 0.0;  /* Minimal density of region points in rectangle. */
  density_th=inputv[2];
  int n_bins = 1024;        /* Number of bins in pseudo-ordering of gradient
                               modulus.                                       */
  return LineSegmentDetection3( n_out, img, X, Y, Z,
                              ang_th, log_eps, density_th, n_bins,
                              NULL,NULL,NULL,NULL ,inputv,inputv_size,inputvorth);
}

/*----------------------------------------------------------------------------*/
/** LSD3 Simple Interface.
 */
double * lsd3b(int * n_out, double * img, int X, int Y, int Z, 
              double *inputv, double inputv_size,double * inputvorth)
{ 
  double ang_th;   /* Gradient angle tolerance in degrees.           */
  ang_th=inputv[4];
  double log_eps = 0.0;     /* Detection threshold: -log10(NFA) > log_eps     */
  log_eps=-log10(inputv[1]);

  double density_th = 0.0;  /* Minimal density of region points in rectangle. */
  density_th=inputv[2];
  int n_bins = 1024;        /* Number of bins in pseudo-ordering of gradient
                               modulus.                                       */
  return LineSegmentDetection3b( n_out, img, X, Y, Z,
                              ang_th, log_eps, density_th, n_bins,
                              NULL,NULL,NULL,NULL ,inputv,inputv_size,inputvorth);
}
double * lsd3bgrad(int * n_out, double * img, int X, int Y, int Z, 
              double *inputv, double inputv_size,double * inputvorth)
{ 
  double ang_th;   /* Gradient angle tolerance in degrees.           */
  ang_th=inputv[4];
  double log_eps = 0.0;     /* Detection threshold: -log10(NFA) > log_eps     */
  log_eps=-log10(inputv[1]);

  double density_th = 0.0;  /* Minimal density of region points in rectangle. */
  density_th=inputv[2];
  int n_bins = 1024;        /* Number of bins in pseudo-ordering of gradient
                               modulus.                                       */
  return LineSegmentDetection3bgrad( n_out, img, X, Y, Z,
                              ang_th, log_eps, density_th, n_bins,
                              NULL,NULL,NULL,NULL ,inputv,inputv_size,inputvorth);
}
/*----------------------------------------------------------------------------*/
/** LSD3 Centerline Simple Interface.
 */
double * lsd3center(int * n_out, 
               double * img, int X, int Y, int Z,
               double * img0, int X0, int Y0, int Z0,
               double *inputv, double inputv_size,double * inputvorth)
{ 
  double ang_th;   /* Gradient angle tolerance in degrees.           */
  ang_th=inputv[4];
  double log_eps = 0.0;     /* Detection threshold: -log10(NFA) > log_eps     */
  log_eps=-log10(inputv[1]);

  double density_th = 0.0;  /* Minimal density of region points in rectangle. */
  density_th=inputv[2];
  int n_bins = 1024;        /* Number of bins in pseudo-ordering of gradient
                               modulus.                                       */
  return LineSegmentDetection3Center( n_out, img, X, Y, Z, img0, X0, Y0, Z0,
                              ang_th, log_eps, density_th, n_bins,
                              NULL,NULL,NULL,NULL ,inputv,inputv_size,inputvorth);
}

/*----------------------------------------------------------------------------*/
/** LSD3 Centerline Simple Interface.
 */
double * lsd3centerb(int * n_out, 
               double * img, int X, int Y, int Z,
               double * img0, int X0, int Y0, int Z0,
               double *inputv, double inputv_size,double * inputvorth,int printreg)
{ 
  double ang_th;   /* Gradient angle tolerance in degrees.           */
  ang_th=inputv[4];
  double log_eps = 0.0;     /* Detection threshold: -log10(NFA) > log_eps     */
  log_eps=-log10(inputv[1]);

  double density_th = 0.0;  /* Minimal density of region points in rectangle. */
  density_th=inputv[2];
  int n_bins = 1024;        /* Number of bins in pseudo-ordering of gradient
                               modulus.                                       */
  return LineSegmentDetection3Centerb( n_out, img, X, Y, Z, img0, X0, Y0, Z0,
                              ang_th, log_eps, density_th, n_bins,
                              NULL,NULL,NULL,NULL ,inputv,inputv_size,inputvorth,printreg);
}

/*----------------------------------------------------------------------------*/
/** LSD+Markov Simple Interface.
 */
double * lsdM(int * n_out, 
               double * img, int X, int Y,
               double * img0, int X0, int Y0,
               double *inputv, double inputv_size)
{
  double ang_th;   /* Gradient angle tolerance in degrees.           */
  ang_th=inputv[4];
  double log_eps = 0.0;     /* Detection threshold: -log10(NFA) > log_eps     */
  log_eps=-log10(inputv[1]);
  
  double density_th = 0.0;  /* Minimal density of region points in rectangle. */
  density_th=inputv[2];
  int n_bins = 1024;        /* Number of bins in pseudo-ordering of gradient
                               modulus.                                       */

  make_markov(img0, X0, Y0, ang_th, n_bins,inputv,inputv_size);                     
  return LineSegmentDetection( n_out, img, X, Y,
                               ang_th, log_eps, density_th, n_bins,
                              NULL,NULL,NULL ,inputv,inputv_size);
}

/*----------------------------------------------------------------------------*/
/** LSD3+Markov Simple Interface.
 */
double * lsd3M(int * n_out, 
               double * img, int X, int Y, int Z,
               double * img0, int X0, int Y0, int Z0,
               double *inputv, double inputv_size,double * inputvorth)
{

  //Instantiate input variables
  double ang_th;   /* Gradient angle tolerance in degrees.           */
  ang_th=inputv[4];
  double log_eps = 0.0;     /* Detection threshold: -log10(NFA) > log_eps     */
  log_eps=-log10(inputv[1]);
  double density_th = 0.0;  /* Minimal density of region points in rectangle. */
  density_th=inputv[2];
  int n_bins = 1024;        /* Number of bins in pseudo-ordering of gradient
                               modulus.                                       */

  //Instantiate constants
  int i,j,il; //iteration variables
    double start,end; //timing variables
  double sizenum = 5.*sqrt(
    (double)X0*(double)X0 + (double)Y0*(double)Y0 + (double)Z0*(double)Z0);
  if(sizenum>pow(10.,6)) sizenum=pow(10.,6);
  inputv[3] = sizenum; //maximum region size
  //printf("IMSIZE: %d\n",(int)sizenum);fflush(stdout);

  //Caluculate "edge feature" Markov kernel from naive image 
  start=omp_get_wtime();
  make_markov3( img0, X0, Y0, Z0, ang_th, n_bins, inputv,inputv_size,0);
  end=omp_get_wtime();
  printf("MAKEMARKOV: %f seconds\n",end-start);fflush(stdout);
  printf("PARAKernel: (p11=%.4f, p10=%.4f), (p11=%.4f, p10=%.4f), (p11=%.4f, p10=%.4f)\n\n",
      inputv[5],inputv[6], inputv[7],inputv[8],inputv[9],inputv[10]);fflush(stdout);
  //t = clock();
  
  //Calculate "centerline feature" Markov kernel from naive image
  start=omp_get_wtime();
  make_markov3( img0, X0, Y0, Z0, ang_th, n_bins, inputvorth,inputv_size,1);
  end=omp_get_wtime();
  printf("MAKEMARKOVORTH: %f seconds\n",end-start);fflush(stdout);
  printf("ORTHKernel: (p11=%.4f, p10=%.4f), (p11=%.4f, p10=%.4f), (p11=%.4f, p10=%.4f)\n\n",
      inputvorth[5],inputvorth[6], inputvorth[7],inputvorth[8],inputvorth[9],inputvorth[10]);fflush(stdout);
  
  //Launch LSD3 algorithm with parallel and orthogonal markov kenrnels

  return LineSegmentDetection3( n_out, img, X, Y, Z,
                                ang_th, log_eps, density_th, n_bins,
                                NULL,NULL,NULL,NULL ,inputv,inputv_size,inputvorth);

}
