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
/*---------------------------------- Import ---------------------------------*/
/*----------------------------------------------------------------------------*/

#ifndef GSOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tgmath.h>
#include <limits.h>
#include <float.h>
#include<string.h>
#include <time.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_qrng.h>
#include<gsl/gsl_sf_gamma.h>
#include<gsl/gsl_sf_trig.h>
#include<sys/mman.h>

#include "rectangles3.h"
#include "constants.h"
#include "misc.h"
#include "tuples.h"
#include "nfa.h"
#include "gradient.h"

/*----------------------------------------------------------------------------*/
/*--------------------------- Rectangle structure ----------------------------*/
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Copy one rectangular prism structure to another.
 */
void rect3_copy(struct rect3 * in, struct rect3 * out)
{
  /* check parameters */
  if( in == NULL || out == NULL ) error("rect3_copy: invalid 'in' or 'out'.");

  /* copy values */
  out->x1 = in->x1;
  out->y1 = in->y1;
  out->z1 = in->z1;
  out->x2 = in->x2;
  out->y2 = in->y2;
  out->z2 = in->z2; 
  out->length = in->length;
  out->width1 = in->width1;
  out->width2 = in->width2;
  out->x = in->x;
  out->y = in->y;
  out->z = in->z;
  out->theta = in->theta;
  for(int i=0;i<3;i++)
  {
    out->dl[i]  = in->dl[i];
    out->dw1[i] = in->dw1[i];
    out->dw2[i] = in->dw2[i];
  }
  out->prec = in->prec;
  out->p = in->p;
  out->density=in->density;
}
/*----------------------------------------------------------------------------*/
/** Rectangle points iterator.

    The integer coordinates of pixels inside a rectangle are
    iteratively explored. This structure keep track of the process and
    functions ri_ini(), ri_inc(), ri_end(), and ri_del() are used in
    the process. An example of how to use the iterator is as follows:
    \code

      struct rect * rec = XXX; // some rectangle
      rect_iter * i;
      for( i=ri_ini(rec); !ri_end(i); ri_inc(i) )
        {
          // your code, using 'i->x' and 'i->y' as coordinates
        }
      ri_del(i); // delete iterator

    \endcode
    The pixels are explored 'column' by 'column', where we call
    'column' a set of pixels with the same x value that are inside the
    rectangle. The following is an schematic representation of a
    rectangle, the 'column' being explored is marked by colons, and
    the current pixel being explored is 'x,y'.

    LSD3 UPDATE: 
    For the prism, ri_ini() is updated to use projections into 
    the polar basis, rather than using ys/ye interpolation.
    We mainain vertex ordering, adjusting so that v4-9 are uniformly higher than 
    v0-3 by the width2 (elevation) distance, so that we iterate over the x/y/z points 
    within the legnth/width/width bounds ('span') of the prism starting from v0, 
    projecting into the nearest coordinate in cartesian space. 
    We maintain the use of rectangular prisms for algorithm extension 
    to thin lines in 3D space, to avoid ground-up reformulations.    
    Hence, detection works on planar features, and data is assumed 
    to have at least one dominant axis in order to establish a 'length'
    concept for the refinement algorithm.   
    \verbatim

                   vx[5],vy[5]vz[5]
                        *   *
                      *       *
                     *           *
                   *               ye
                  *                :  *
             vx[4],vy[4]vz[4]      :     *
                    *              :        *
                       *          x,y,z        *
                          *        :              *
                             *     :            vx[6],vy[6]vz[6]
                                *  :                *
                                   ys              *
                                  (ze)*           *
                                   :    *       *
                                  :        *   *
              vx[1],vy[1]vz[1]   :    vx[7],vy[7],vz[7]
                 *   *          :
                *       *      :
               *           *  (zs)
              *               ye
             *                :  *
        vx[0],vy[0]vz[0]      :     *
               *              :        *
            xs    *          x,y,z        *
                :    *        :              *
                   :    *     :            vx[2],vy[2]vz[2]
                      :   *   :                *
        y               :     ys              *
        ^                  :     *           *
        |                    :      *       *
        |                      xe      *   *
        +---> x                   vx[3],vy[3]vz[3]
       -
      -
     z

    \endverbatim
    The first 'column' to be explored is the one with the smaller x
    value. Each 'column' is explored starting from the pixel of the
    'column' (inside the rectangle) with the smallest y value.

    The four corners of the rectangle are stored in order that rotates
    around the corners at the arrays 'vx[]' and 'vy[]'. The first
    point is always the one with smaller x value.

    'x' and 'y' are the coordinates of the pixel being explored. 'ys'
    and 'ye' are the start and end values of the current column being
    explored. So, 'ys' < 'ye'.
 */



/*----------------------------------------------------------------------------*/
/** Free memory used by a 3D rectangle iterator.
 */
void ri3_del(rect3_iter * iter)
{
  if( iter == NULL ) error("ri_del: NULL iterator.");
  free( (void *) iter );
}

/*----------------------------------------------------------------------------*/
/** Check if the 3D iterator finished the full iteration.
 */
int ri3_end(rect3_iter * i)
{
  /* check input */
  if( i == NULL ) error("ri3_end: NULL iterator.");
  /* reached end of projected space*/
  return (i->zt > i->zspan) && (i->xt > i->xspan) && (i->yt > i->yspan);
}

/*----------------------------------------------------------------------------*/
/** Project the coordinates from spherical to cartsian space, with origin at v0.
 *  Forming the rotation matrix from the spherical bases, Q=[n,t] 
 *  and given an input vector vt incrementing along the polar bases in the polar frame, 
 *  we recover the pixel coordinates in the cartesian frame by : 
 *  (vd) = v0 + Q(vt) 
 *   */
void up_all3(rect3_iter * i)
{
 /*
 double xshift = ((double)(i->xt)*(i->dl[0])) 
    + ((double)(i->yt)*(i->dl[1])) + ((double)(i->zt)*(i->dl[2]));

 double yshift = ((double)(i->xt)*(i->dw1[0])) 
    + ((double)(i->yt)*(i->dw1[1])) + ((double)(i->zt)*(i->dw1[2]));

  double zshift = ((double)(i->xt)*(i->dw2[0])) 
    + ((double)(i->yt)*(i->dw2[1])) + ((double)(i->zt)*(i->dw2[2]));
  */
 double xshift = ((double)(i->xt)*(i->dl[0])) 
    + ((double)(i->yt)*(i->dw1[0])) + ((double)(i->zt)*(i->dw2[0]));

 double yshift = ((double)(i->xt)*(i->dl[1])) 
    + ((double)(i->yt)*(i->dw1[1])) + ((double)(i->zt)*(i->dw2[1]));

  double zshift = ((double)(i->xt)*(i->dl[2])) 
    + ((double)(i->yt)*(i->dw1[2])) + ((double)(i->zt)*(i->dw2[2]));

 // int xtemp=(int)ceil(vx[0]+((double)i->xspan * dl[0] + (double)i->yspan * dw1[0] + (double)i->zspan * dw2[0])-1);


 int xmax = (int) ceil(i->vx[0]+((double)(i->xspan)*(i->dl[0])) 
    + ((double)(i->yspan)*(i->dw1[0])) + ((double)(i->zspan)*(i->dw2[0])))-1;

 int ymax = (int)ceil(i->vy[0]+((double)(i->xspan)*(i->dl[1])) 
    + ((double)(i->yspan)*(i->dw1[1])) + ((double)(i->zspan)*(i->dw2[1])));

  int zmax = (int) ceil(i->vz[0]+((double)(i->xspan)*(i->dl[2])) 
    + ((double)(i->yspan)*(i->dw1[2])) + ((double)(i->zspan)*(i->dw2[2])));
  
  //x projection 
  i->xd = i->vx[0]+xshift;
   if (i->x != (int) ceil(i->xd)-1) {
    i->x = (int) ceil(i->xd)-1;
    i->update = 1;}
  //y projection 
  i->yd = i->vy[0] +yshift;
  if (i->y != (int) ceil(i->yd)) {
    i->y = (int) ceil(i->yd);
    i->update = 1;}
  //z projection 
  i->zd = i->vz[0] +zshift;
  if (i->z != (int) ceil(i->zd)) {
    i->z = (int) ceil(i->zd);
    i->update = 1;}
  if(1==0){//ri3_end(i)){
  //printf("LWW %d, %d, %d; MAX;  SHIFT %.2f, %.2f, %.2f\n",i->xspan,i->yspan,i->zspan,xshift,yshift,zshift);fflush(stdout);
  //printf("v0:  [%.2f, %.2f, %.2f]; v4:[%.2f, %.2f, %.2f]; vmax: [%d, %d, %d]; vend[%d, %d, %d]\n",
  //i->vx[0],i->vy[0],i->vz[0],
  //i->vx[7],i->vy[7],i->vz[7],
  //xmax,ymax,zmax,
  //i->x, i->y, i->z
  //);fflush(stdout);
  //printf("v0:  [%d, %d, %d]; v4:[%d, %d, %d]; vmax: [%d, %d, %d]; vend[%d, %d, %d]\n",
  //(int)ceil(i->vx[0])-1,(int)ceil(i->vy[0]),(int)ceil(i->vz[0]),
  //(int)ceil(i->vx[7])-1,(int)ceil(i->vy[7]),(int)ceil(i->vz[7]),
  //xmax,ymax,zmax,
  //i->x, i->y, i->z
  //);fflush(stdout);
//if ((int)ceil(i->vz[7]>100)) error("BIG BOY");
//if (i->vz[0]>i->vz[7]) error("DIR ERR");
}
}

/*----------------------------------------------------------------------------*/
/** Iterate pixels in the polar reference frame for integer step sizes, 
 *  projecting into cartesian coordinates from the origin v0.*/
void ri3_inc(rect3_iter * i)
{
  /* check input */
  if( i == NULL ) error("ri_inc: NULL iterator.");
  i->update=0;
  //Increment integer-wise in the rotated plane 
  while(i->update==0)
  {
    //3D: New OUTER layer
    // if not at end of exploration, inc. z to next 'depth'
    if(!(ri3_end(i))) i->zt++;
    //if at end of col, inc. x to next row
    while( (i->zt > i->zspan) && !(ri3_end(i)) )
    {
      //if not at end of exploration, inc. y to next 'col'
      if(!(ri3_end(i))) i->yt++;
      while( (i->yt > i->yspan) && !(ri3_end(i)) )   
      {
        i->xt++;
        if(ri3_end(i)){up_all3(i); return;}
        i->yt=0;
      }
      if(ri3_end(i)) {up_all3(i); return;}
      i->zt = 0;
    }
    //Update cartesian pixel coordinate 
    up_all3(i);
  }
}

/*----------------------------------------------------------------------------*/
/** Create and initialize a rectangle iterator in 3D. 
 * */
rect3_iter * ri3_ini(struct rect3 * r)
{
  double vx[8],vy[8],vz[8];
  int n,offset;
  rect3_iter * i;

  /* check parameters */
  if( r == NULL ) error("ri3_ini: invalid rectangle.");

  /* get memory */
  i = (rect3_iter *) malloc(sizeof(rect3_iter));
  if( i == NULL ) error("ri3_ini: Not enough memory.");

  /* build list of rectangle corners ordered
     in a circular way around the rectangle 
     with two z level sets (v0-v3 on lower z)
     NOTE: width1 is along the azimithal direction, 
           width2 is along the elevation direction,
  */
  //printf("\nR: (%.2f, %.2f, %.2f ) - (%.2f, %.2f, %.2f) for LWW(%.2f, %.2f, %.2f) \n",r->x1,r->y1,r->z1, r->x2, r->y2, r->z2,r->length,r->width1,r->width2);fflush(stdout);
  //printf("Jacobian: [%.2f, %.2f, %.2f], [%.2f, %.2f, %.2f], [%.2f, %.2f, %.2f]\n",
  //r->dl[0],r->dl[1], r->dl[2],  r->dw1[0],r->dw1[1], r->dw1[2],  r->dw2[0],r->dw2[1], r->dw2[2]
//);fflush(stdout);
  vx[0] = r->x1 - (r->dw1[0] * r->width1 / 2.0) - (r->dw2[0] * r->width2 / 2.0);
  vy[0] = r->y1 - (r->dw1[1] * r->width1 / 2.0) - (r->dw2[1] * r->width2 / 2.0);
  vz[0] = r->z1 - (r->dw1[2] * r->width1 / 2.0) - (r->dw2[2] * r->width2 / 2.0);

  vx[1] = r->x2 - (r->dw1[0] * r->width1 / 2.0) - (r->dw2[0] * r->width2 / 2.0);
  vy[1] = r->y2 - (r->dw1[1] * r->width1 / 2.0) - (r->dw2[1] * r->width2 / 2.0);
  vz[1] = r->z2 - (r->dw1[2] * r->width1 / 2.0) - (r->dw2[2] * r->width2 / 2.0);

  vx[2] = r->x2 + (r->dw1[0] * r->width1 / 2.0) - (r->dw2[0] * r->width2 / 2.0);
  vy[2] = r->y2 + (r->dw1[1] * r->width1 / 2.0) - (r->dw2[1] * r->width2 / 2.0);
  vz[2] = r->z2 + (r->dw1[2] * r->width1 / 2.0) - (r->dw2[2] * r->width2 / 2.0);

  vx[3] = r->x1 + (r->dw1[0] * r->width1 / 2.0) - (r->dw2[0] * r->width2 / 2.0);
  vy[3] = r->y1 + (r->dw1[1] * r->width1 / 2.0) - (r->dw2[1] * r->width2 / 2.0);
  vz[3] = r->z1 + (r->dw1[2] * r->width1 / 2.0) - (r->dw2[2] * r->width2 / 2.0);

  //upper z
  vx[4] = r->x1 - (r->dw1[0] * r->width1 / 2.0) + (r->dw2[0] * r->width2 / 2.0);
  vy[4] = r->y1 - (r->dw1[1] * r->width1 / 2.0) + (r->dw2[1] * r->width2 / 2.0);
  vz[4] = r->z1 - (r->dw1[2] * r->width1 / 2.0) + (r->dw2[2] * r->width2 / 2.0);

  vx[5] = r->x2 - (r->dw1[0] * r->width1 / 2.0) + (r->dw2[0] * r->width2 / 2.0);
  vy[5] = r->y2 - (r->dw1[1] * r->width1 / 2.0) + (r->dw2[1] * r->width2 / 2.0);
  vz[5] = r->z2 - (r->dw1[2] * r->width1 / 2.0) + (r->dw2[2] * r->width2 / 2.0);

  vx[6] = r->x2 + (r->dw1[0] * r->width1 / 2.0) + (r->dw2[0] * r->width2 / 2.0);
  vy[6] = r->y2 + (r->dw1[1] * r->width1 / 2.0) + (r->dw2[1] * r->width2 / 2.0);
  vz[6] = r->z2 + (r->dw1[2] * r->width1 / 2.0) + (r->dw2[2] * r->width2 / 2.0);

  vx[7] = r->x1 + (r->dw1[0] * r->width1 / 2.0) + (r->dw2[0] * r->width2 / 2.0);
  vy[7] = r->y1 + (r->dw1[1] * r->width1 / 2.0) + (r->dw2[1] * r->width2 / 2.0);
  vz[7] = r->z1 + (r->dw1[2] * r->width1 / 2.0) + (r->dw2[2] * r->width2 / 2.0);
  
  /* compute rotation of index of corners needed so that the first
     point has the smaller x from among the smaller z.
     if one side is vertical, thus two corners have the same smaller z
     value, the one with the largest x value is selected as the first 
     if one side is vertical, thus two corners have the same smaller x
     value, the one with the largest y value is selected as the first.
   
     UPDATE: cycle over dw1 plane and then dw2 plane in order to cycle 
             through to the proper v0 relations
   */
  /* 
  offset=1;
  double tx[8],ty[8], tz[8];
  while(vx[0]>vx[1] || vx[0]>vx[2] || vx[0]>vx[3])
  {
    while(vx[0]>vx[1] || vx[0]>vx[2] || vx[0]>vx[3])
    {
      for(n=0;n<8;n++)
      {
        tx[n]=vx[n];
        ty[n]=vy[n];
        tz[n]=vz[n];
      }
      for(n=0;n<4;n++)
      {
        //z1
        vx[n] = tx[(offset+n)%4];
        vy[n] = ty[(offset+n)%4];
        vz[n] = tz[(offset+n)%4];
        //z2
        vx[n+4] = tx[(offset+n)%4+4];
        vy[n+4] = ty[(offset+n)%4+4];
        vz[n+4] = tz[(offset+n)%4+4];
      }
    }
    if (vx[0]>vx[4])
    {
      for(n=0;n<8;n++)
      {
        tx[n]=vx[n];
        ty[n]=vy[n];
        tz[n]=vz[n];
      }
      for(n=0;n<4;n++)
      {
        //z1
        vx[n] = tx[n+4];
        vy[n] = ty[n+4];
        vz[n] = tz[n+4];
        //z2
        vx[n+4] = tx[n];
        vy[n+4] = ty[n];
        vz[n+4] = tz[n];
      }
    }
  }
  */
  for(n=0;n<8;n++)
  {
      i->vx[n]=vx[n];
      i->vy[n]=vy[n];
      i->vz[n]=vz[n];
  }
  /* Set an initial condition.

     The values are set to values that will cause 'ri_inc' (that will
     be called immediately) to initialize correctly the first 'column'
     and compute the limits 'ys' and 'ye'.

     'y' is set to the integer value of vy[0], the starting corner.

     'ys' and 'ye' are set to very small values, so 'ri_inc' will
     notice that it needs to start a new 'column'.

     The smallest integer coordinate inside of the rectangle is
     'ceil(vx[0])'. The current 'x' value is set to that value minus
     one, so 'ri_inc' (that will increase x by one) will advance to
     the first 'column'.
   */
  i->x = (int) ceil(i->vx[0])-1; //innermost index, analogous to 2D
  i->y = (int) ceil(i->vy[0]);
  i->z = (int) ceil(i->vz[0]);
  i->xspan = r->length; 
  i->yspan = r->width1;
  i->zspan = r->width2;
  i->xt = 0; i->yt = 0; i->zt = 0; i->xd = 0; i->yd = 0; i->zd = 0;
  /*Set absolute-angles in order to ensure proper x-incrementation */
  double saz,caz,sel,cel;
  sincos(r->theta->az,&saz,&caz);
  sincos(r->theta->el,&sel,&cel);
  double  dl[3] =  {caz*sel, saz*sel,  cel};
  double dw1[3] =  {-saz, caz, 0.};
  double dw2[3] =  {caz*cel, saz*cel, -sel};  
  



 double xmax = i->vx[0]+((double)(i->xspan)*(i->dl[0])) 
    + ((double)(i->yspan)*(i->dw1[0])) + ((double)(i->zspan)*(i->dw2[0]));

 double ymax = i->vy[0]+((double)(i->xspan)*(i->dl[1])) 
    + ((double)(i->yspan)*(i->dw1[1])) + ((double)(i->zspan)*(i->dw2[1]));

 double zmax = i->vz[0]+((double)(i->xspan)*(i->dl[2])) 
    + ((double)(i->yspan)*(i->dw1[2])) + ((double)(i->zspan)*(i->dw2[2]));
  





  
  //saz=fabs(saz);caz=fabs(caz);sel=fabs(sel);cel=fabs(cel); 
  double tempaz,tempel;
  double tsaz,tcaz,tsel,tcel;
  /*
  if(((double)i->xspan * dl[0] + (double)i->yspan * dl[1] + (double)i->zspan * dl[2]) <0.)
  {
    tempaz = r->theta->az;
    tempel = r->theta->el;
    align3(&tempaz,&tempel);
    sincos(tempaz,&tsaz,&tcaz);
    sincos(tempel,&tsel,&tcel);
    dl[0] = tcaz*tsel;
    dl[1] = tsaz*tsel;
    dl[2] = tcel;
   }
  */
  
  if(i->vy[0] > i->vy[3]) 
  //if((((double)i->xspan * dw1[0] + (double)i->yspan * dw1[1]+ ( double)i->zspan * dw1[2]) <0.) && (vy[0]<vy[2]) )
  {
    tempaz = r->theta->az;
    tempel = r->theta->el;
    align3(&tempaz,&tempel);
    sincos(tempaz,&tsaz,&tcaz);
    sincos(tempel,&tsel,&tcel);
    //dw1[0] =  -tsaz;
    //dw1[1] =  tcaz;
    //dw1[0] *=-1.;
    //dw1[1] *=-1.;

  } 
  if(i->vz[0] > i->vz[4]) 
  //if((((double)i->xspan * dw2[0] + (double)i->yspan * dw2[1] + (double)i->zspan * dw2[2]) <0.) &&(vz[0]<vz[4]) )
  {
    tempaz = r->theta->az;
    tempel = r->theta->el;
    align3(&tempaz,&tempel);
    sincos(tempaz,&tsaz,&tcaz);
    sincos(tempel,&tsel,&tcel);
    //dw2[0] =  tcaz*tcel;  
   // dw2[1] =  tsaz*tcel;  
    //dw2[2] =  -tsel;  
    //dw2[0]*=-1.;
    //dw2[1]*=-1.;
    //dw2[2]*=-1.;
  }

  int xxx=0;
  for(int ix=0; ix<8;ix++)
  {
    if (vx[ix]>vx[xxx]) xxx=ix;
  }  
  //printf("\n idx %d",xxx);fflush(stdout);
 
  //if(((double)i->xspan * dl[0] + (double)i->yspan * dl[1] + (double)i->zspan * dl[2]) <0.)
  double dx = ((double)i->xspan * dl[0] + (double)i->yspan * dw1[0] + (double)i->zspan * dw2[0]);
  //if (vx[0]<vx[xxx])
  if (1==0)//abs(xtemp-vx[4])<2
  {
    tempaz = r->theta->az;
    tempel = r->theta->el;
    align3(&tempaz,&tempel);
    sincos(tempaz,&tsaz,&tcaz);
    sincos(tempel,&tsel,&tcel);
    //dl[0] = tcaz*tsel;
    //dl[1] = tsaz*tsel;
    //dl[2] = tcel;
    //dw1[0] =  -tsaz;
    //dw1[1] =  tcaz;
    //dw2[0] =  tcaz*tcel;  
    //dw2[1] =  tsaz*tcel;  
    //dw2[2] =  -tsel;  
    dl[0] *=1.;
    dl[1] *=1.;
    dl[2] *=1.;
    dw1[0] *=-1;
    dw1[1] *=-1.;
    dw2[0] *=-1.;  
    dw2[1] *=-1;  
    dw2[2] *=-1.;  
   }

 /*
  //if(i->vy[0] > i->vy[2]) 
  //if((((double)i->xspan * dw1[0] + (double)i->yspan * dw1[1]+ ( double)i->zspan * dw1[2]) <0.) && (vy[0]<vy[3]) )
  if((((double)i->xspan * dl[1] + (double)i->yspan * dw1[1]+ ( double)i->zspan * dw2[1]) <0.) && (vy[0]<vy[3]) )
  {
    tempaz = r->theta->az;
    tempel = r->theta->el;
    align3(&tempaz,&tempel);
    sincos(tempaz,&tsaz,&tcaz);
    sincos(tempel,&tsel,&tcel);
    //dw1[0] =  -tsaz;
    //dw1[1] =  tcaz;
    dw1[0] =  tcaz*tcel;  
    dw1[1] =  tsaz*tcel;  
    dw1[2] =  -tsel;  

  } 
  //if(i->vz[0] > i->vz[5]) 
  //if((((double)i->xspan * dw2[0] + (double)i->yspan * dw2[1] + (double)i->zspan * dw2[2]) <0.) &&(vz[0]<vz[5]) )
  if((((double)i->xspan * dl[2] + (double)i->yspan * dw1[2] + (double)i->zspan * dw2[2]) <0.) &&(vz[0]<vz[5]) )
  {
    tempaz = r->theta->az;
    tempel = r->theta->el;
    align3(&tempaz,&tempel);
    sincos(tempaz,&tsaz,&tcaz);
    sincos(tempel,&tsel,&tcel);
    dw2[0] =  -tsaz;
    dw2[1] =  tcaz;
    //dw2[0] =  tcaz*tcel;  
    //dw2[1] =  tsaz*tcel;  
   */ 
    
    
    
    
  /* Instantiate unit normal and tangent vectors, 
   * fixing the normal vector to the 1st quadrent 
   * to ensure proper exploration along the x (length) vector
   * from the origin at v0, the smallest x coordinate
   */

  for(n=0;n<3;n++) 
  {
    i->dl[n]  = dl[n];
    i->dw1[n] = dw1[n];
    i->dw2[n] = dw2[n];
  }
  /* Correct the orientation of the tangent axis as needed
   * to ensure valid exploration, since v0 should have 
   * the smallest y value, and v2 the largest.  So to for Z 
   * */
  //if(i->vy[0] < i->vy[2]) for(n=0;n<3;n++) i->dw1[n]*=-1.; 
  //if(i->vz[0] < i->vz[5]) for(n=0;n<3;n++) i->dw2[n]*=-1.;
  /* advance to the first pixel */
  ri3_inc(i);
  return i;
}

/*----------------------------------------------------------------------------*/
/** Compute a rectangle's NFA for the rectangular prism.
 *
 *  Ideally, we would grow points using parallel alignment of gradient vectors,
 *  and then find the NFA by counting alignments orthogonal to the priciple axis.
 *  This would encourage detection through the principle axis of a line, rather than
 *  individual detections on each surface face. 
 *
 *  To avoid statistical complexities from moving between parallel/orthogonal frameworks,
 *  we choose instead to count alignments parallel to each principl axis individually, 
 *  and then determine the NFA from the maximum alignment count. 
 *  This equates the statistics between region growing and NFA calculation, 
 *  and accounts for alignments across a wide 2D surface as opposed to a 1D line edge.
 *
 *  UPDATE: 
 *    Test counting along the intermediate axis only.
 *  TODO: 
 *    Apply new orthogonal mnfa for true-orthogonal consideration?    
 *    Note that we are growing in the parallel framework, and switching to orthogonal kernel may be detrimental 
 *
 */
double rect3_nfa(struct rect3 * rec, grads angles, double logNT,double *image,int N,int minreg)
{
  rect3_iter * i;
  int pts = 0;
  int alg = 0;

  /* check parameters */
  if( rec == NULL ) error("rect3_nfa: invalid rectangle.");
  if( angles == NULL ) error("rect3_nfa: invalid grads structure 'angles'.");
  if( angles->az->data ==NULL || angles->el->data ==NULL) error("rect3_nfa:invalid az/el images withing grads structure 'angles'");

  int tester=0;
  double grads_az, grads_el;
  int xsize=(int) angles->az->xsize;
  int ysize=(int) angles->az->ysize;
  int zsize=(int) angles->az->zsize;
  /* local alias to principal components */
  double theta_az = rec->theta->az;
  double theta_el = rec->theta->el;
  double theta_az2 = rec->theta->az2;
  double theta_el2 = rec->theta->el2;
  double theta_az3 = rec->theta->az3;
  double theta_el3 = rec->theta->el3;
  double prec = rec->prec;
  int alg2=0;
  int alg3=0;

  if (prec<0.0) error("rect3_nfa: 'prec' must be positive");

  /* precompute trig functions for principal axes */
  double cprec = cos(prec);
  double sGel,cGel,sTel,cTel,sTel2,cTel2,sTel3,cTel3;
  sincos(theta_el,&sTel,&cTel);
  sincos(theta_el2,&sTel2,&cTel2);
  sincos(theta_el3,&sTel3,&cTel3);
  double diff,diff2,diff3;

  /* compute the total number of pixels and of aligned points in 'rec' */
  for(i=ri3_ini(rec); !ri3_end(i); ri3_inc(i)) // rectangle iterator 
  {
    if( i->x >= 0 && i->y >= 0 && i->z >=0 &&
          i->x < xsize && i->y < ysize && i->z < zsize)
    {
      ++pts;
      grads_az = angles->az->data[ i->z + zsize*(i->x + i->y * xsize) ]; 
      grads_el = angles->el->data[ i->z + zsize*(i->x + i->y * xsize) ];
      
      //in-line version of 'isaligned3' to reduce function calls
      if( !(grads_az == NOTDEF) && !(grads_el == NOTDEF) ) 
      {
        if( (isfinite(grads_az)) && (isfinite(grads_el)) ) 
        {
          if(isaligned3ORTH(grads_az,grads_el,theta_az,theta_el,prec)) 
            ++alg;
        }
      }    

    }
  }

  ri3_del(i); /* delete iterator */
  rec->density=(double)alg/(double)pts;
  if (logNT==0) return -1;
  double nfaout=0;
  if(pts<minreg)
      nfaout= -1;
  else
    nfaout= nfa(pts,alg,rec->p,logNT,image,N); /* compute NFA value */
  //WARNING IN CASE OF ERROR: 
  if (isinf(nfaout))
  {
    double imgval =  image[alg*N+pts];
    printf("\tWARNING: rect3_nfa: non-finite nfa: k/pts = %d / %d, mnfa = %.2e\n",alg,pts,imgval);fflush(stdout);
  }
  return nfaout;
}

/*----------------------------------------------------------------------------*/
/** Compute a rectangle's NFA for the rectangular prism.
 *
 *  For centerline detection, we have a prior centerline estimate and orthogonal propability kernels
 *  TODO: 
 *    Merge with rect3_nfa and test effects.  
 *    rect3_nfa lacks *pset, but we can reintroduce those.  Transfer should be direct and easy, but 
 *    need to confirm efficacy since GROWTH is parallel, and may invalidate binomial law.
 *
 */
double rect3_nfaORTH(struct rect3 * rec, grads angles, double logNT,double *image,double *pset, int N,int minreg)
{
  rect3_iter * i;
  int pts = 0;
  int alg = 0;

  /* check parameters */
  if( rec == NULL ) error("rect3_nfa: invalid rectangle.");
  if( angles == NULL ) error("rect3_nfa: invalid grads structure 'angles'.");
  if( angles->az->data ==NULL || angles->el->data ==NULL) error("rect3_nfa:invalid az/el images withing grads structure 'angles'");

  int tester=0;
  double grads_az, grads_el;
  int xsize=(int) angles->az->xsize;
  int ysize=(int) angles->az->ysize;
  int zsize=(int) angles->az->zsize;
  /* move major, intermediate, and minor axis orientations to local memory */
  double theta_az = rec->theta->az;
  double theta_el = rec->theta->el;
  double theta_az2 = rec->theta->az2;
  double theta_el2 = rec->theta->el2;
  double theta_az3 = rec->theta->az3;
  double theta_el3 = rec->theta->el3;
  double prec = rec->prec;
  int alg2=0;
  int alg3=0;
  if (prec<0.0) error("rect3_nfa: 'prec' must be positive");

  /* precompute trig functions for principal axes */
  double cprec = cos(prec);
  double sGel,cGel,sTel,cTel,sTel2,cTel2,sTel3,cTel3;
  sincos(theta_el,&sTel,&cTel);
  sincos(theta_el2,&sTel2,&cTel2);
  sincos(theta_el3,&sTel3,&cTel3);
  double diff,diff2,diff3;
  /* compute the total number of pixels and of aligned points in 'rec' */
  for(i=ri3_ini(rec); !ri3_end(i); ri3_inc(i)) // rectangle iterator 
  {
    if( i->x >= 0 && i->y >= 0 && i->z >=0 &&
          i->x < xsize && i->y < ysize && i->z < zsize)
    {
      ++pts;
      grads_az = angles->az->data[ i->z + zsize*(i->x + i->y * xsize) ]; 
      grads_el = angles->el->data[ i->z + zsize*(i->x + i->y * xsize) ];
      
      //in-line version of 'isaligned3' to reduce function calls
      if( !(grads_az == NOTDEF) && !(grads_el == NOTDEF) ) 
      {
        if( (isfinite(grads_az)) && (isfinite(grads_el)) ) 
        {
          if(isaligned3ORTH(grads_az,grads_el,theta_az,theta_el,prec)) 
            ++alg;
        }
      }    
     

     }
  }
  ri3_del(i); /* delete iterator */
  rec->density=(double)alg/(double)pts;
  if (logNT==0) return -1;
  //save nfa for switching to binomial estimation if outside markov table or nonfinite
  double nfaout=0.;
  double tempnfaout=0.;
  double fixnfaout=0.;
  double divout=0;
  if(pts<minreg)
    nfaout=-1;
  else
  {
    if (pts<N)
    {
      nfaout=nfa(pts,alg,rec->p,logNT,image,N); /* compute NFA value */
      //tempnfaout=nfaORTH(pts,alg,rec->p,logNT,pset,N,&divout); /* compute NFA value */
      if(!isfinite(nfaout))
        {nfaout=nfaORTH(pts,alg,rec->p,logNT,pset,N,&divout); /* compute NFA value */
         //printf("ERROR: RECT3: INF NFA "); fflush(stdout);
        }
    }
    else
      {
      nfaout = nfaORTH(pts,alg,rec->p,logNT,pset,N,&divout); /* compute NFA value */
   
      }
  }
 return nfaout;
}
