/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/constants.h: constants definition

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
/*---------------------------- constants.h --------------------------------*/
/*----------------------------------------------------------------------------*/

//Open header
#ifndef CONSTANTS_HEADER
#define CONSTANTS_HEADER

/*--------------------------------------------------------------------------
  Define constants 
*/
/** ln(10) */
#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif /* !M_LN10 */

/** PI */
#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif /* !M_PI */

#ifndef RADIANS_TO_DEGREES
#define RADIANS_TO_DEGREES (180.0/M_PI)
#endif /*!RADIANS_TO_DEGREES*/

#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

/** Label for pixels with undefined gradient. */
#ifndef NOTDEF
#define NOTDEF -1024.0
#endif
/** 3/2 pi */
#ifndef M_3_2_PI
#define M_3_2_PI 4.71238898038
#endif
/** 2 pi */
#ifndef M_2__PI
#define M_2__PI  6.28318530718
#endif
/** Label for pixels not used in yet. */
#ifndef NOTUSED
#define NOTUSED 0
#endif
/** Label for pixels already used in detection. */
#ifndef USED
#define USED    1
#endif 

/*----------------------------------------------------------------------------*/
/** Doubles relative error factor
 */
#ifndef RELATIVE_ERROR_FACTOR
#define RELATIVE_ERROR_FACTOR 100.0
#endif

//Close header
#endif /* !LSD_HEADER */