/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/regions2.h - region approximation, measurement, and improvement functions for the 2D case.

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
/*---------------------------- regions2.h --------------------------------*/
/*----------------------------------------------------------------------------*/

//Open header
#ifndef REGIONS2_HEADER
#define REGIONS2_HEADER

#include "tuples.h"
#include "rectangles2.h"
struct point; //incomplete forward declaration from misc;


//Define functions
double get_theta( struct point * reg, int reg_size, double x, double y,
                        image_double modgrad, double reg_angle, double prec );
void region2rect( struct point * reg, int reg_size,
                        image_double modgrad, double reg_angle,
                        double prec, double p, struct rect * rec );
void region2rectORTH( struct point * reg, int reg_size,
                        image_double modgrad, double reg_angle,
                        double prec, double p, struct rect * rec );
void region_grow( int x, int y, image_double angles, struct point * reg,
                        int * reg_size, double * reg_angle, image_char used,
                        double prec );
void region_growORTH( int x, int y, image_double angles, 
                        struct point * reg, int * reg_size, 
                        double * reg_angle, double * lstheta, 
                        image_char used, double prec);

double rect_improve( struct rect * rec, image_double angles,
                        double logNT, double log_eps,
                        double* mnfa,double* mnfa_2,double* mnfa_4,
                        int Nnfa,int minsize, int minsize2,int minsize4 );
double rect_improve_update(struct rect  r, image_double angles,double logNT,int Nnfa,
                        double* mnfa, double* mnfap, int minsize,
                        double* mnfa_2,double* mnfap_2, int minsize2,
                        double* mnfa_4,double* mnfap_4, int minsize4,
                        double p1check, double p2check,
                        struct rect * rec,double log_nfa,int orth);
double rect_improveORTH( struct rect * rec, image_double angles,
                        double logNT, double log_eps,
                        double* mnfa,double* mnfa_2,double* mnfa_4,
                        double*mnfap,double*mnfap_2,double*mnfap_4,
                        int Nnfa,int minsize, int minsize2,int minsize4, int orth );
int reduce_region_radius( struct point * reg, int * reg_size,
                        image_double modgrad, double reg_angle,
                        double prec, double p, struct rect * rec,
                        image_char used, image_double angles,
                        double density_th );
int reduce_region_radiusORTH( struct point * reg, int * reg_size,
                        image_double modgrad, double reg_angle,
                        double prec, double p, struct rect * rec,
                        image_char used, image_double angles,
                        double density_th );
int refine( struct point * reg, int * reg_size, image_double modgrad,
                        double reg_angle, double prec, double p, struct rect * rec,
                        image_char used, image_double angles, double density_th );
int refineORTH( struct point * reg, int * reg_size, image_double modgrad,
                   double reg_angle, double prec, double p, struct rect * rec,
                   image_char used, image_double angles, double density_th ,int orth);







//Close header
#endif
