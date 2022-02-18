/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/rectangles3.c - iterator and nfa functions for the 3D case.


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
/*---------------------------- rectangles3.h --------------------------------*/
/*----------------------------------------------------------------------------*/

//Open header
#ifndef RECTANGLES3_HEADER
#define RECTANGLES3_HEADER

#include "misc.h" //required for anlges
#include "gradient.h"

//Define functions
/*----------------------------------------------------------------------------*/
/** Rectangular prism structure: line segment with two orthogonal widths.
 */
struct rect3
{
  double x1,y1,z1,x2,y2,z2;     /* first and second point of the line segment */
  double length,width1,width2;  /* rectangle width */
  double x,y,z;                 /* center of the rectangle */
  angles3 theta;                /* az/el angle as struct angle3 */
  double dl[3],dw1[3],dw2[3];   /* spherical basis vectors of the principal axis*/
  double prec;                  /* tolerance angle */
  double p;                     /* probability of a point with angle within 'prec' */
  double density;
};

void rect3_copy(struct rect3 * in, struct rect3 * out);

/*See rectangles3.c*/
typedef struct 
{
  double vx[8];  /* rectangle's corner X coordinates in circular order */
  double vy[8];  /* rectangle's corner Y coordinates in circular order */
  double vz[8]; 
  int x,y,z; // pixel coordinates in original image frame
  //INTRODUCED FOR LSD3, for projected iteration in polar bases 
  int update;              // indicator for if projected pixel is a new coordinate
  double xd,yd,zd;         // projected coordinate in cartesian basis
  int xt,yt,zt;            // explored coordinate in sphereical basis
  int xspan, yspan, zspan; // range of explorable space in spherical basis
  double dl[3];            //vector for rotating the x coordinate (normal vector)
  double dw1[3];           //vector for rotating the y coordinate (azimuth tangent)
  double dw2[3];           //vector for rotating the z coordinate (elevation tangent)
} rect3_iter;

void ri3_del(rect3_iter * iter);
int ri3_end(rect3_iter * i);
void up_all3(rect3_iter * i);
void ri3_inc(rect3_iter * i);
rect3_iter * ri3_ini(struct rect3 * r);
double rect3_nfa(struct rect3 * rec, grads angles, 
                        double logNT,double *image,int N,int minreg);
double rect3_nfaORTH(struct rect3 * rec, grads angles, 
                            double logNT,double *image,double *pset, int N,int minreg);


//Close header
#endif
