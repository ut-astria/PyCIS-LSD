/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/markov.h: Determine first-order markov kernel for angular tolerance options

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
#
------------------------------------------------------------------------------*/   

/*----------------------------------------------------------------------------*/
/*---------------------------- markov.h --------------------------------*/
/*----------------------------------------------------------------------------*/

//Open header
#ifndef MARKOV_HEADER
#define MARKOV_HEADER

//Define functions
void NFA_matrix(double *output,double p0,double p11,double p01,int N);
void make_markov( double * img, int X, int Y,
                           double ang_th, int n_bins,
                           double * inputv,double inputv_size);
void make_markovORTH( double * img, int X, int Y,
                           double ang_th, int n_bins,
                           double * inputv,double inputv_size);
static int isaligned3_markovVORTH(double grads_az,double grads_el,double cprec);
static int isaligned3_markovHORTH(double grads_az,double grads_el,double cprec);
static int isaligned3_markovDORTH(double grads_az,double grads_el,double cprec);
static int isaligned3_markovV(double grads_az,double grads_el,double cprec);
static int isaligned3_markovH(double grads_az,double grads_el,double cprec);
static int isaligned3_markovD(double grads_az,double grads_el,double cprec);
void make_markov3( double * img, int X, int Y, int Z,
                          double ang_th, int n_bins, double * inputv,double inputv_size, 
                          int orth);
                          
//Close header
#endif
