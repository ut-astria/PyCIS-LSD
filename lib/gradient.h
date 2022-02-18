/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/gradient.h: compute gradient image using Sobel or Gradient-by-Ratio in 2D/3D

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
/*---------------------------- gradient.h --------------------------------*/
/*----------------------------------------------------------------------------*/

//Open header
#ifndef GRADIENT_HEADER
#define GRADIENT_HEADER

#include "tuples.h"
struct coorlist; //incomplete forward declaration 
struct coorlist3; //incomplete forward declaration 


//Define Functions
image_double ll_angle( image_double in,
                              struct coorlist ** list_p, void ** mem_p,
                              image_double * modgrad, unsigned int n_bins,double alpha);

/*----------------------------------------------------------------------------*/
/* Storage structure for 3D gradients. grads->az, grads->el,
 * instatiate as struct grads newgrad*/
typedef struct grads_s {
  image3_double az;
  image3_double el;
} * grads;

grads new_grads(unsigned int xsize, unsigned int ysize, unsigned int zsize);
void free_grads(grads i);
grads ll_angle3( image3_double in,
                        struct coorlist3 ** list_p, void ** mem_p,
                        image3_double * modgrad, 
                        unsigned int n_bins,double alpha);

int isaligned( int x, int y, image_double angles, double theta,  double prec );
int isalignedORTH( int x, int y, image_double angles, double theta,  double prec );
int isaligned3(double grads_az,double grads_el,double theta_az,double theta_el,double prec);
int isaligned3_sign(double grads_az,double grads_el,double theta_az,double theta_el,double prec);
int isaligned3ORTH(double grads_az,double grads_el,double theta_az,double theta_el,double prec);
void align3(double * az, double * el);
//Close header
#endif /* !LSD_HEADER */