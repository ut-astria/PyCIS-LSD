/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/misc.h: helper functions

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
/*---------------------------- misc.h --------------------------------*/
/*----------------------------------------------------------------------------*/

//Open header
#ifndef MISC_HEADER
#define MISC_HEADER

//Define functions
int max1(int x, int y);
int min1(int x, int y);
double max2(double x, double y);
double min2(double x, double y);
/*----------------------------------------------------------------------------*/
/** Chained list of coordinates.
 */
struct coorlist
{
  int x,y;
  struct coorlist * next;
};
struct coorlist3
{
    int x,y,z;
    struct coorlist3 * next;
};
/*----------------------------------------------------------------------------*/
/** A point (or pixel).
 */
struct point {int x,y;};
struct point3 {int x,y,z;};

void error(char * msg);
int double_equal(double a, double b);

double dist(double x1, double y1, double x2, double y2);
double dist3(double x1, double y1, double z1, double x2, double y2, double z2);
double line_angle(double x1, double y1, double x2, double y2);

/*----------------------------------------------------------------------------*/
/*Orientation of a line in sphereical coordinates
 * (az,el) represents the principal axis of a line, while
 * (az3,el3) and (az2,el2) represent the minor and intermediate axes, 
 * for consideration of planar events.
*/
typedef struct angles3_s
{
  double az,el;  
  double az2,el2;
  double az3,el3;
} * angles3;

angles3 new_angles3(double az, double el);
void free_angles3(angles3 i);
angles3 line_angle3(double x1, double y1, double z1, 
        double x2, double y2, double z2);
double angle_diff(double a, double b);
double angle_diff_signed(double a, double b);

//Close header
#endif
