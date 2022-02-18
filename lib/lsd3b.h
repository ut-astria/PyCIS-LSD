/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/lsd3b.h: main pipline for 1st order line detection 
  2D/3D integration for temporal aliasing
  edge line / center line two-layer detection for narrow 3D linear features

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
/*---------------------------- lsd3b.h --------------------------------*/
/*----------------------------------------------------------------------------*/

//Open header
#ifndef LSD3B_HEADER
#define LSD3B_HEADER

#include "tuples.h"

double * LineSegmentDetectionb( int * n_out,
                               double * img, int X, int Y,
                               double ang_th, double log_eps, double density_th,
                               int n_bins, int ** reg_img, 
                               int * reg_x, int * reg_y ,double * inputv, double inputv_size,
                               unsigned char ** used_img,
                               double * output, double * output_2, double * output_4,double logNT);

double * LineSegmentDetectionCenterb( int * n_out,
                               double * img, int X, int Y,
                               double ang_th, double log_eps, double density_th,
                               int n_bins, int ** reg_img, 
                               int * reg_x, int * reg_y ,double * inputv, double inputv_size,
                               unsigned char ** used_img,
                               double * output, double * output_2, double * output_4,
                               double * outputP, double * outputP_2, double * outputP_4,
                               ntuple_list out, int k,double logNT);

double * LineSegmentDetection3b( int * n_out,
                               double * img, int X, int Y, int Z,
                               double ang_th, double log_eps, double density_th,
                               int n_bins, int ** reg_img, 
                               int * reg_x, int * reg_y, int * reg_z, 
                               double * inputv, double inputv_size, double * inputvorth);

double * LineSegmentDetection3bgrad( int * n_out,
                               double * img, int X, int Y, int Z,
                               double ang_th, double log_eps, double density_th,
                               int n_bins, int ** reg_img, 
                               int * reg_x, int * reg_y, int * reg_z, 
                               double * inputv, double inputv_size, double * inputvorth);

double * LineSegmentDetection3Centerb( int * n_out,
                               double * img, int X, int Y, int Z,
                               double * img0, int X0, int Y0, int Z0,
                               double ang_th, double log_eps, double density_th,
                               int n_bins, int ** reg_img, 
                               int * reg_x, int * reg_y, int * reg_z, 
                               double * inputv, double inputv_size, double * inputvorth,
				int printreg);

#endif
