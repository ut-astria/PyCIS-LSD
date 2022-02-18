/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/regions3.h - region approximation, measurement, and improvement functions for the 3D case.

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
/*---------------------------- regions3.h --------------------------------*/
/*----------------------------------------------------------------------------*/

//Open header
#ifndef REGIONS3_HEADER
#define REGIONS3_HEADER

#include "misc.h" //required for angles3
#include "tuples.h"
#include "gradient.h" //regquired fro grads
#include "rectangles3.h"
struct point3; //incomplete forward declaration from misc.h

//Define functions
angles3 get_theta3( struct point3 * reg, int reg_size, double x, double y, double z,
                        image3_double modgrad, angles3 reg_angle, double prec, int orth );
void region2rect3( struct point3 * reg, int reg_size,
                        image3_double modgrad, angles3 reg_angle,
                        double prec, double p, struct rect3 * rec , int orth );
void region3_grow(int x, int y,int z, grads angles, 
                        struct point3 * reg,
                        int * reg_size, angles3 * reg_angle, 
                        image3_char used,double prec ,int NOUT);
void region3_growORTH(int x, int y,int z, 
                        image3_double modgrad, grads angles, 
                        struct point3 * reg, int * reg_size, 
                        angles3 * reg_angle,  angles3 * lstheta, 
                            image3_char used,double prec ,int NOUT);
double rect3_improve_update(struct rect3  r, grads angles,double logNT,int Nnfa,
                        double* mnfa, double* mnfap, int minsize,
                        double* mnfa_2,double* mnfap_2, int minsize2,
                        double* mnfa_4,double* mnfap_4, int minsize4,
                        double p1check, double p2check,
                        struct rect3 * rec,double log_nfa,int orth);
double rect3_improve( struct rect3 * rec, grads angles,
                        double logNT, double log_eps,
                        double* mnfa,double* mnfa_2,double* mnfa_4,
                        double*mnfap,double*mnfap_2,double*mnfap_4,
                        int Nnfa,int minsize, int minsize2,int minsize4,int orth);
int reduce_region3_radius( struct point3 * reg, int * reg_size,
                        image3_double modgrad, angles3 reg_angle,
                        double prec, double p, struct rect3 * rec,
                        image3_char used, grads angles,
                        double density_th , int orth);
int refine3( struct point3 * reg, int * reg_size, image3_double modgrad,
                        angles3 reg_angle, double prec, double p, struct rect3 * rec,
                        image3_char used, grads angles,
                        double density_th , int NOUT, int orth);

                                

//Close header
#endif
