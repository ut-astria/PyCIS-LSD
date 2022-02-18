/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/rectangles2.c - iterator and nfa functions for the 2D case.

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
/*---------------------------- rectangles2.h --------------------------------*/
/*----------------------------------------------------------------------------*/

//Open header
#ifndef RECTANGLES2_HEADER
#define RECTANGLES2_HEADER

#include "tuples.h"

//Define functions
/*----------------------------------------------------------------------------*/
/** Rectangle structure: line segment with width.
 */
struct rect
{
  double x1,y1,x2,y2;  /* first and second point of the line segment */
  double width;        /* rectangle width */
  double x,y;          /* center of the rectangle */
  double theta;        /* angle */
  double dx,dy;        /* (dx,dy) is vector oriented as the line segment */
  double prec;         /* tolerance angle */
  double p;            /* probability of a point with angle within 'prec' */
};

void rect_copy(struct rect * in, struct rect * out);

//see rectangles2.c
typedef struct 
{
double vx[4];     // rectangle's corner X coordinates in circular order 
double vy[4];     // rectangle's corner Y coordinates in circular order 
int x,y;          // coordinates of currently explored pixel 
double ys,ye;     //LSDSAR: endpoins of y column at x step
//INTRODUCED FOR LSD3, for projected iteration in polar bases 
int update;       // indicator for if projected pixel is a new coordinate
double xd,yd;     // projected coordinate in cartesian basis
int xt,yt;        // explored coordinate in sphereical basis
int xspan, yspan; // range of explorable space in spherical basis
double dl[2];     //vector for rotating the x coordinate (normal vector)
double dn[2];     //vector for rotating the y coordinate (tangent vector)
} rect_iter;

void ri_del(rect_iter * iter);
int ri_end(rect_iter * i);
void up_all(rect_iter * i);
void ri_inc(rect_iter * i);
rect_iter * ri_ini(struct rect * r);
double rect_nfa(struct rect * rec, image_double angles, 
                        double logNT,double *image,int N,int minreg);
double rect_nfaORTH(struct rect * rec, image_double angles, 
                        double logNT,double *image,double *pset,int N,int minreg);


//Close header
#endif
