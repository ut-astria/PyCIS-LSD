/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/misc.c: helper functions

Benjamin Feuge-Miller: benjamin.g.miller@utexas.edu
The University of Texas at Austin, 
Oden Institute Computational Astronautical Sciences and Technologies (CAST) group
*Date of Modification: September 03, 2021

#
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
#

------------------------------------------------------------------------------*/   

/*----------------------------------------------------------------------------*/
/*---------------------------------- Import ---------------------------------*/
/*----------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tgmath.h>
#include <limits.h>
#include <float.h>
#include<string.h>
#include <time.h>
#include <gsl/gsl_sf_trig.h>
#include <sys/mman.h>

#include "misc.h"
#include "constants.h"


/*----------------------------------------------------------------------------*/
/*------------------------- Miscellaneous functions --------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Max/min helper functions
 */
int max1(int x, int y)
{ return x<y?y:x;}
int min1(int x, int y)
{return x<y?x:y;}

double max2(double x, double y)
{ return x<y?y:x;}
double min2(double x, double y)
{return x<y?x:y;}


/*----------------------------------------------------------------------------*/
/** Fatal error, print a message to standard-error output and exit.
 */
void error(char * msg)
{
  fprintf(stderr,"LSD Error: %s\n",msg);
  exit(EXIT_FAILURE);
}



/*----------------------------------------------------------------------------*/
/** Compare doubles by relative error.

    The resulting rounding error after floating point computations
    depend on the specific operations done. The same number computed by
    different algorithms could present different rounding errors. For a
    useful comparison, an estimation of the relative rounding error
    should be considered and compared to a factor times EPS. The factor
    should be related to the cumulated rounding error in the chain of
    computation. Here, as a simplification, a fixed factor is used.
 */
int double_equal(double a, double b)
{
  double abs_diff,aa,bb,abs_max;

  /* trivial case */
  if( a == b ) return TRUE;

  abs_diff = fabs(a-b);
  aa = fabs(a);
  bb = fabs(b);
  abs_max = aa > bb ? aa : bb;

  /* DBL_MIN is the smallest normalized number, thus, the smallest
     number whose relative error is bounded by DBL_EPSILON. For
     smaller numbers, the same quantization steps as for DBL_MIN
     are used. Then, for smaller numbers, a meaningful "relative"
     error should be computed by dividing the difference by DBL_MIN. */
  if( abs_max < DBL_MIN ) abs_max = DBL_MIN;

  /* equal if relative error <= factor x eps */
  return (abs_diff / abs_max) <= (RELATIVE_ERROR_FACTOR * DBL_EPSILON);
}


/*----------------------------------------------------------------------------*/
/*---------------------------- Line functions --------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Computes Euclidean distance between point (x1,y1) and point (x2,y2).
 */
double dist(double x1, double y1, double x2, double y2)
{
  return sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
}

/*----------------------------------------------------------------------------*/
/** Computes Euclidean distance between point (x1,y1,z1) and point (x2,y2,z2).
 */
double dist3(double x1, double y1, double z1, double x2, double y2, double z2)
{
  return sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) );
}

/*----------------------------------------------------------------------------*/
/* Orientation of a line*/
double line_angle(double x1, double y1, double x2, double y2)
{
  return atan2( y2-y1 , x2-x1 );
}
/*----------------------------------------------------------------------------*/
/** Instantiate a new angles object given a principal axis orientation.*/
angles3 new_angles3(double az, double el)
{
  angles3 image;
  /* get memory */
  image = (angles3) malloc( sizeof(struct angles3_s) );
  if( image == NULL ) error("not enough memory.");
  //input
  image->az = az;
  image->el = el; 
  return image;
}


/*----------------------------------------------------------------------------*/
/** Free memory used in angles3 */
void free_angles3(angles3 i)
{
  if( i == NULL)
    error("free_angles3: invalid input image.");
  free( (void *) i );
}

/*----------------------------------------------------------------------------*/
/*Given a unit Cartesian line, create and angles structure of the polar orientation */
angles3 line_angle3(double x1, double y1, double z1, 
        double x2, double y2, double z2)
{
    double az =  atan2( y2-y1 , x2-x1) ;
    double el =  acos((z2-z1)/dist3(x1,y1,z1,x2,y2,z2));
    angles3 azel=new_angles3(az,el);
    return azel;
}
/*----------------------------------------------------------------------------*/
/** Absolute value angle difference.
 */
double angle_diff(double a, double b)
{
  a -= b;
  while( a <= -M_PI ) a += M_2__PI;
  while( a >   M_PI ) a -= M_2__PI;
  if( a < 0.0 ) a = -a;
  return a;
}

/*----------------------------------------------------------------------------*/
/** Signed angle difference.
 */
double angle_diff_signed(double a, double b)
{
  a -= b;
  while( a <= -M_PI ) a += M_2__PI;
  while( a >   M_PI ) a -= M_2__PI;
  return a;
}

