/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/regions2.c - region approximation, measurement, and improvement functions for the 2D case.

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



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tgmath.h>
#include <limits.h>
#include <float.h>
#include<string.h>
#include <time.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_sf_gamma.h>
#include<gsl/gsl_sf_trig.h>
#include<sys/mman.h>

#include "regions2.h"
#include "constants.h"
#include "misc.h"
#include "tuples.h"
#include "nfa.h"
#include "rectangles2.h"
#include "gradient.h"

/*----------------------------------------------------------------------------*/
/*---------------------------------- Regions ---------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Compute region's angle as the principal inertia axis of the region.

    The following is the region inertia matrix A:
    @f[
        A = \left(\begin{array}{cc}
              Ixx & Ixy \\
              Ixy & Iyy \\
             \end{array}\right)
    @f]
    where
      Ixx =   sum_i G(i).(y_i - cx)^2
      Iyy =   sum_i G(i).(x_i - cy)^2
      Ixy = - sum_i G(i).(x_i - cx).(y_i - cy)
    and
    - G(i) is the gradient norm at pixel i, used as pixel's weight.
    - x_i and y_i are the coordinates of pixel i.
    - cx and cy are the coordinates of the center of th region.

    lambda1 and lambda2 are the eigenvalues of matrix A,
    with lambda1 >= lambda2. They are found by solving the
    characteristic polynomial:

      det( lambda I - A) = 0

    that gives:

      lambda1 = ( Ixx + Iyy + sqrt( (Ixx-Iyy)^2 + 4.0*Ixy*Ixy) ) / 2
      lambda2 = ( Ixx + Iyy - sqrt( (Ixx-Iyy)^2 + 4.0*Ixy*Ixy) ) / 2

    To get the line segment direction we want to get the angle the
    eigenvector associated to the smallest eigenvalue. We have
    to solve for a,b in:

      a.Ixx + b.Ixy = a.lambda2
      a.Ixy + b.Iyy = b.lambda2

    We want the angle theta = atan(b/a). It can be computed with
    any of the two equations:

      theta = atan( (lambda2-Ixx) / Ixy )
    or
      theta = atan( Ixy / (lambda2-Iyy) )

    When |Ixx| > |Iyy| we use the first, otherwise the second (just to
    get better numeric precision).
 */
double get_theta( struct point * reg, int reg_size, double x, double y,
                         image_double modgrad, double reg_angle, double prec )
{
  double lambda,theta,weight;
  double Ixx = 0.0;
  double Iyy = 0.0;
  double Ixy = 0.0;
  int i;

  /* check parameters */
  if( reg == NULL ) error("get_theta: invalid region.");
  if( reg_size <= 1 ) error("get_theta: region size <= 1.");
  if( modgrad == NULL || modgrad->data == NULL )
    error("get_theta: invalid 'modgrad'.");
  if( prec < 0.0 ) error("get_theta: 'prec' must be positive.");

  /* compute inertia matrix */
  for(i=0; i<reg_size; i++)
    {
      weight = modgrad->data[ reg[i].x + reg[i].y * modgrad->xsize ];
      Ixx += ( (double) reg[i].y - y ) * ( (double) reg[i].y - y ) * weight;
      Iyy += ( (double) reg[i].x - x ) * ( (double) reg[i].x - x ) * weight;
      Ixy -= ( (double) reg[i].x - x ) * ( (double) reg[i].y - y ) * weight;
    }
  if( double_equal(Ixx,0.0) && double_equal(Iyy,0.0) && double_equal(Ixy,0.0) )
    error("get_theta: null inertia matrix.");

  /* compute smallest eigenvalue */
  lambda = 0.5 * ( Ixx + Iyy - sqrt( (Ixx-Iyy)*(Ixx-Iyy) + 4.0*Ixy*Ixy ) );

  /* compute angle */
  theta = fabs(Ixx)>fabs(Iyy) ? atan2(lambda-Ixx,Ixy) : atan2(Ixy,lambda-Iyy);

  /* The previous procedure doesn't cares about orientation,
     so it could be wrong by 180 degrees. Here is corrected if necessary. */
  if( angle_diff(theta,reg_angle) > prec ) theta += M_PI;

  return theta;
}

/*----------------------------------------------------------------------------*/
/** Computes a rectangle that covers a region of points.
 */
void region2rect( struct point * reg, int reg_size,
                         image_double modgrad, double reg_angle,
                         double prec, double p, struct rect * rec )
{
  double x,y,dx,dy,l,w,theta,weight,sum,l_min,l_max,w_min,w_max;
  int i;

  /* check parameters */
  if( reg == NULL ) error("region2rect: invalid region.");
  if( reg_size <= 1 ) error("region2rect: region size <= 1.");
  if( modgrad == NULL || modgrad->data == NULL )
    error("region2rect: invalid image 'modgrad'.");
  if( rec == NULL ) error("region2rect: invalid 'rec'.");

  /* center of the region:
     It is computed as the weighted sum of the coordinates
     of all the pixels in the region. The norm of the gradient
     is used as the weight of a pixel. The sum is as follows:
       cx = \sum_i G(i).x_i
       cy = \sum_i G(i).y_i
     where G(i) is the norm of the gradient of pixel i
     and x_i,y_i are its coordinates.
   */
  x = y = sum = 0.0;
  for(i=0; i<reg_size; i++)
    {
      weight = modgrad->data[ reg[i].x + reg[i].y * modgrad->xsize ];
      x += (double) reg[i].x * weight;
      y += (double) reg[i].y * weight;
      sum += weight;
    }
  if( sum <= 0.0 ) error("region2rect: weights sum equal to zero.");
  x /= sum;
  y /= sum;

  /* theta */
  theta = get_theta(reg,reg_size,x,y,modgrad,reg_angle,prec);

  /* length and width:

     'l' and 'w' are computed as the distance from the center of the
     region to pixel i, projected along the rectangle axis (dx,dy) and
     to the orthogonal axis (-dy,dx), respectively.

     The length of the rectangle goes from l_min to l_max, where l_min
     and l_max are the minimum and maximum values of l in the region.
     Analogously, the width is selected from w_min to w_max, where
     w_min and w_max are the minimum and maximum of w for the pixels
     in the region.
   */
  dx = cos(theta);
  dy = sin(theta);
  l_min = l_max = w_min = w_max = 0.0;
  for(i=0; i<reg_size; i++)
    {
      l =  ( (double) reg[i].x - x) * dx + ( (double) reg[i].y - y) * dy;
      w = -( (double) reg[i].x - x) * dy + ( (double) reg[i].y - y) * dx;

      if( l > l_max ) l_max = l;
      if( l < l_min ) l_min = l;
      if( w > w_max ) w_max = w;
      if( w < w_min ) w_min = w;
    }

  /* store values */
  rec->x1 = x + l_min * dx;
  rec->y1 = y + l_min * dy;
  rec->x2 = x + l_max * dx;
  rec->y2 = y + l_max * dy;
  rec->width = w_max - w_min;
  rec->x = x;
  rec->y = y;
  rec->theta = theta;
  rec->dx = dx;
  rec->dy = dy;
  rec->prec = prec;
  rec->p = p;

  /* we impose a minimal width of one pixel

     A sharp horizontal or vertical step would produce a perfectly
     horizontal or vertical region. The width computed would be
     zero. But that corresponds to a one pixels width transition in
     the image.
   */
  if( rec->width < 1.0 ) rec->width = 1.0;
}
void region2rectORTH( struct point * reg, int reg_size,
                         image_double modgrad, double reg_angle,
                         double prec, double p, struct rect * rec )
{
  double x,y,dx,dy,l,w,theta,weight,sum,l_min,l_max,w_min,w_max;
  int i;

  /* check parameters */
  if( reg == NULL ) error("region2rect: invalid region.");
  if( reg_size <= 1 ) error("region2rect: region size <= 1.");
  if( modgrad == NULL || modgrad->data == NULL )
    error("region2rect: invalid image 'modgrad'.");
  if( rec == NULL ) error("region2rect: invalid 'rec'.");

  /* center of the region:
     It is computed as the weighted sum of the coordinates
     of all the pixels in the region. The norm of the gradient
     is used as the weight of a pixel. The sum is as follows:
       cx = \sum_i G(i).x_i
       cy = \sum_i G(i).y_i
     where G(i) is the norm of the gradient of pixel i
     and x_i,y_i are its coordinates.
   */
  x = y = sum = 0.0;
  for(i=0; i<reg_size; i++)
    {
      weight = modgrad->data[ reg[i].x + reg[i].y * modgrad->xsize ];
      x += (double) reg[i].x * weight;
      y += (double) reg[i].y * weight;
      sum += weight;
    }
  if( sum <= 0.0 ) error("region2rect: weights sum equal to zero.");
  x /= sum;
  y /= sum;

  /* theta */
  theta = reg_angle;//get_theta(reg,reg_size,x,y,modgrad,reg_angle,prec);

  /* length and width:

     'l' and 'w' are computed as the distance from the center of the
     region to pixel i, projected along the rectangle axis (dx,dy) and
     to the orthogonal axis (-dy,dx), respectively.

     The length of the rectangle goes from l_min to l_max, where l_min
     and l_max are the minimum and maximum values of l in the region.
     Analogously, the width is selected from w_min to w_max, where
     w_min and w_max are the minimum and maximum of w for the pixels
     in the region.
   */
  dx = cos(theta);
  dy = sin(theta);
  l_min = l_max = w_min = w_max = 0.0;
  for(i=0; i<reg_size; i++)
    {
      l =  ( (double) reg[i].x - x) * dx + ( (double) reg[i].y - y) * dy;
      w = -( (double) reg[i].x - x) * dy + ( (double) reg[i].y - y) * dx;

      if( l > l_max ) l_max = l;
      if( l < l_min ) l_min = l;
      if( w > w_max ) w_max = w;
      if( w < w_min ) w_min = w;
    }

  /* store values */
  rec->x1 = x + l_min * dx;
  rec->y1 = y + l_min * dy;
  rec->x2 = x + l_max * dx;
  rec->y2 = y + l_max * dy;
  rec->width = w_max - w_min;
  rec->x = x;
  rec->y = y;
  rec->theta = theta;
  rec->dx = dx;
  rec->dy = dy;
  rec->prec = prec;
  rec->p = p;

  /* we impose a minimal width of one pixel

     A sharp horizontal or vertical step would produce a perfectly
     horizontal or vertical region. The width computed would be
     zero. But that corresponds to a one pixels width transition in
     the image.
   */
  if( rec->width < 1.0 ) rec->width = 1.0;
}

/*----------------------------------------------------------------------------*/
/** Build a region of pixels that share the same angle, up to a
    tolerance 'prec', starting at point (x,y).
 */
void region_grow( int x, int y, image_double angles, struct point * reg,
                         int * reg_size, double * reg_angle, image_char used,
                         double prec )
{
  double sumdx,sumdy;
  int xx,yy,i;

  /* check parameters */
  if( x < 0 || y < 0 || x >= (int) angles->xsize || y >= (int) angles->ysize )
    error("region_grow: (x,y) out of the image.");
  if( angles == NULL || angles->data == NULL )
    error("region_grow: invalid image 'angles'.");
  if( reg == NULL ) error("region_grow: invalid 'reg'.");
  if( reg_size == NULL ) error("region_grow: invalid pointer 'reg_size'.");
  if( reg_angle == NULL ) error("region_grow: invalid pointer 'reg_angle'.");
  if( used == NULL || used->data == NULL )
    error("region_grow: invalid image 'used'.");

  /* first point of the region */
  *reg_size = 1;
  reg[0].x = x;
  reg[0].y = y;
  *reg_angle = angles->data[x+y*angles->xsize];  /* region's angle */
  sumdx = cos(*reg_angle);
  sumdy = sin(*reg_angle);
  used->data[x+y*used->xsize] = USED;

  /* try neighbors as new region points */
  for(i=0; i<*reg_size; i++)
    for(xx=reg[i].x-1; xx<=reg[i].x+1; xx++)
      for(yy=reg[i].y-1; yy<=reg[i].y+1; yy++)
        if( xx>=0 && yy>=0 && xx<(int)used->xsize && yy<(int)used->ysize &&
            used->data[xx+yy*used->xsize] != USED &&
            isaligned(xx,yy,angles,*reg_angle,prec) )
          {
            //printf("grown...");fflush(stdout);
            /* add point */
            used->data[xx+yy*used->xsize] = USED;
            reg[*reg_size].x = xx;
            reg[*reg_size].y = yy;
            ++(*reg_size);

            /* update region's angle */
            sumdx += cos( angles->data[xx+yy*angles->xsize] );
            sumdy += sin( angles->data[xx+yy*angles->xsize] );
            *reg_angle = atan2(sumdy,sumdx);
          }
}

/*----------------------------------------------------------------------------*/
/** Build a region of pixels that share the same angle,
 * orthogonal to a known prior centerline, 
 *  up to a tolerance 'prec', starting at point (x,y).
 * NOTE: This is much easier than the 3D case, since
 * the orthogonal directions are finite in the spatial plane
 */
void region_growORTH( int x, int y, image_double angles, 
                        struct point * reg, int * reg_size, 
                        double * reg_angle, double * lstheta, 
                        image_char used, double prec)
{
  double sumdx,sumdy;
  int xx,yy,i;

  /* check parameters */
  if( x < 0 || y < 0 || x >= (int) angles->xsize || y >= (int) angles->ysize )
    error("region_grow: (x,y) out of the image.");
  if( angles == NULL || angles->data == NULL )
    error("region_grow: invalid image 'angles'.");
  if( reg == NULL ) error("region_grow: invalid 'reg'.");
  if( reg_size == NULL ) error("region_grow: invalid pointer 'reg_size'.");
  if( reg_angle == NULL ) error("region_grow: invalid pointer 'reg_angle'.");
  if( used == NULL || used->data == NULL )
    error("region_grow: invalid image 'used'.");

  /* first point of the region */
  *reg_size = 1;
  reg[0].x = x;
  reg[0].y = y;


  //*reg_angle = angles->data[x+y*angles->xsize];  /* region's angle */
  //Define prior as the centerline, so sumdx/dy should be parallel
  *reg_angle = *lstheta;
  sumdx = cos(*reg_angle);
  sumdy = sin(*reg_angle);
  used->data[x+y*used->xsize] = USED;
  /*To check 'orthogonal alignment', can simply add pi/2V*/
  double  temp_angle;

  /* try neighbors as new region points */
  for(i=0; i<(*reg_size); i++)
    for(xx=reg[i].x-1; xx<=reg[i].x+1; xx++)
      for(yy=reg[i].y-1; yy<=reg[i].y+1; yy++)
      {

        if( xx>=0 && yy>=0 && xx<(int)used->xsize && yy<(int)used->ysize &&
            used->data[xx+yy*used->xsize] != USED &&
            isalignedORTH(xx,yy,angles,*reg_angle,prec) )
          {
            /* add point */
            used->data[xx+yy*used->xsize] = USED;
            reg[*reg_size].x = xx;
            reg[*reg_size].y = yy;
            ++(*reg_size);

            /* update region's centerline (orthogonal) angle */
            /*ensure orthogonal direction is properly aligned parallel to region*/
            temp_angle = angles->data[xx+yy*angles->xsize];
            if ((sumdx*(cos(temp_angle)) + sumdy*sin(temp_angle))<0)
            {
              temp_angle = temp_angle + M_PI/2.;
              if (temp_angle>M_PI) temp_angle -= M_PI;
            }
            sumdx += cos( temp_angle );
            sumdy += sin( temp_angle );
            //fit prior line in region growing - we know these have a high probability of existance
            //*reg_angle = atan2(sumdy,sumdx);
          }
      }
}

/*----------------------------------------------------------------------------*/
/** Try some rectangles variations to improve NFA value. Only if the
    rectangle is not meaningful (i.e., log_nfa <= log_eps).
 */
double rect_improve( struct rect * rec, image_double angles,
                            double logNT, double log_eps,
                            double* mnfa,double* mnfa_2,double* mnfa_4,
                            int Nnfa,int minsize, int minsize2,int minsize4 )
{
  struct rect r;
  double log_nfa,log_nfa_new;
  double delta = 0.5;
  double delta_2 = delta / 2.0;
  double factor = 1./sqrt(2.);
  int n;
  rect_copy(rec,&r);
  if(r.p>0.1)
    log_nfa = rect_nfa(&r,angles,logNT,mnfa,Nnfa,minsize);
  else if(r.p>0.05)
     log_nfa= rect_nfa(&r,angles,logNT,mnfa_2,Nnfa,minsize2);
  else
     log_nfa = rect_nfa(&r,angles,logNT,mnfa_4,Nnfa,minsize4);
  
  if( log_nfa > log_eps ) return log_nfa;

  /* try finer precisions */
  rect_copy(rec,&r);
  for(n=0; n<1; n++)
  {
    r.p /= 2.0;
    r.prec = r.p * M_PI;
    if(r.p>0.1)
      log_nfa_new = rect_nfa(&r,angles,logNT,mnfa,Nnfa,minsize);
    else if(r.p>0.05)
         log_nfa_new = rect_nfa(&r,angles,logNT,mnfa_2,Nnfa,minsize2);
    else
         log_nfa_new = rect_nfa(&r,angles,logNT,mnfa_4,Nnfa,minsize4);
    if( log_nfa_new > log_nfa )
    {
      log_nfa = log_nfa_new;
      rect_copy(&r,rec);
    }
  }
  if( log_nfa > log_eps ) return log_nfa;

  /* try to reduce width */
  rect_copy(rec,&r);
  for(n=0; n<5; n++)
  {
    //if( (r.width - delta) >= 0.5 )
    if((r.width*factor) >= 0.5) //( (r.width - delta) >= 0.5 )
    {
      //r.width -= delta;
      r.width *= factor;
      if(r.p>0.1)
        log_nfa_new = rect_nfa(&r,angles,logNT,mnfa,Nnfa,minsize);
      else if(r.p>0.05)
        log_nfa_new = rect_nfa(&r,angles,logNT,mnfa_2,Nnfa,minsize2);
      else
        log_nfa_new = rect_nfa(&r,angles,logNT,mnfa_4,Nnfa,minsize4);
      if( log_nfa_new > log_nfa )
      {
        rect_copy(&r,rec);
        log_nfa = log_nfa_new;
      }
    }
  }
  if( log_nfa > log_eps ) return log_nfa;

  /* try to reduce one side of the rectangle */
  rect_copy(rec,&r);
  for(n=0; n<5; n++)
  {
    //if( (r.width - delta) >= 0.5 )
    //{
    if((r.width*factor) >= 0.5) //( (r.width - delta) >= 0.5 )
    {
      delta_2 = r.width*(1.0-factor)/2.;
      r.x1 += -r.dy * delta_2;
      r.y1 +=  r.dx * delta_2;
      r.x2 += -r.dy * delta_2;
      r.y2 +=  r.dx * delta_2;
      //r.width -= delta;
      r.width *= factor;
       
      if(r.p>0.1)
        log_nfa_new = rect_nfa(&r,angles,logNT,mnfa,Nnfa,minsize);
      else if(r.p>0.05)
        log_nfa_new = rect_nfa(&r,angles,logNT,mnfa_2,Nnfa,minsize2);
      else
         log_nfa_new = rect_nfa(&r,angles,logNT,mnfa_4,Nnfa,minsize4);
       
      if( log_nfa_new > log_nfa )
      {
        rect_copy(&r,rec);
        log_nfa = log_nfa_new;
      }
    }
  }
  if( log_nfa > log_eps ) return log_nfa;

  /* try to reduce the other side of the rectangle */
  rect_copy(rec,&r);
  for(n=0; n<5; n++)
  {
    //if( (r.width - delta) >= 0.5 )
    //{
    if((r.width*factor) >= 0.5) //( (r.width - delta) >= 0.5 )
    {
      delta_2 = r.width*(1.0-factor)/2.;
      r.x1 -= -r.dy * delta_2;
      r.y1 -=  r.dx * delta_2;
      r.x2 -= -r.dy * delta_2;
      r.y2 -=  r.dx * delta_2;
      //r.width -= delta;
      r.width *= factor;
      if(r.p>0.1)
        log_nfa_new = rect_nfa(&r,angles,logNT,mnfa,Nnfa,minsize);
      else if(r.p>0.05)
        log_nfa_new = rect_nfa(&r,angles,logNT,mnfa_2,Nnfa,minsize2);
      else
        log_nfa_new = rect_nfa(&r,angles,logNT,mnfa_4,Nnfa,minsize4);
      if( log_nfa_new > log_nfa )
      {
        rect_copy(&r,rec);
        log_nfa = log_nfa_new;
      }
    }
  }
  if( log_nfa > log_eps ) return log_nfa;

  /* try even finer precisions */
  rect_copy(rec,&r);
  for(n=0; n<1; n++)
  {
    r.p /= 2.0;
    r.prec = r.p * M_PI;
    if(r.p>0.1)
      log_nfa_new = rect_nfa(&r,angles,logNT,mnfa,Nnfa,minsize);
    else if(r.p>0.05)
      log_nfa_new = rect_nfa(&r,angles,logNT,mnfa_2,Nnfa,minsize2);
    else
      log_nfa_new = rect_nfa(&r,angles,logNT,mnfa_4,Nnfa,minsize4);
    if( log_nfa_new > log_nfa )
    {
      log_nfa = log_nfa_new;
      rect_copy(&r,rec);
    }
  }
  return log_nfa;
}

/*----------------------------------------------------------------------------*/
/** Helper function for trying variations to improve NFA value.
 */
double rect_improve_update(struct rect  r, image_double angles,double logNT,int Nnfa,
                                 double* mnfa, double* mnfap, int minsize,
                                 double* mnfa_2,double* mnfap_2, int minsize2,
                                 double* mnfa_4,double* mnfap_4, int minsize4,
                                 double p1check, double p2check,
                                 struct rect * rec,double log_nfa,int orth)
{
  double log_nfa_new;
  if (dist(rec->x1,rec->y1,rec->x2,rec->y2)<=rec->width) return -101.;
  /*Pick which NFA function to call based on sequence flag "orth" and the r.p value*/
  if (orth==0)
  {
    if(r.p>=p1check)
      log_nfa_new = rect_nfa(&r,angles,logNT,mnfa,Nnfa,minsize);
    else if(r.p>=p2check)
      log_nfa_new = rect_nfa(&r,angles,logNT,mnfa_2,Nnfa,minsize2);
    else
      log_nfa_new = rect_nfa(&r,angles,logNT,mnfa_4,Nnfa,minsize4);
  }
  else if (orth==1)
  {
    if(r.p>=p1check)
      log_nfa_new = rect_nfaORTH(&r,angles,logNT,mnfa,mnfap,Nnfa,minsize);
    else if(r.p>=p2check)
      log_nfa_new = rect_nfaORTH(&r,angles,logNT,mnfa_2,mnfap_2,Nnfa,minsize2);
    else
      log_nfa_new = rect_nfaORTH(&r,angles,logNT,mnfa_4,mnfap_4,Nnfa,minsize4);
  }
  else error("rect_improve_update: post flag not recognized 0/1");

  /*Return original estimate or update rectangle appropriately*/
  if (log_nfa==-1)
    return log_nfa_new;
  if(log_nfa_new > log_nfa)
  {
    rect_copy(&r,rec);
    return log_nfa_new;
  }
  else return log_nfa;
   
}

/*----------------------------------------------------------------------------*/
/** Try some rectangles variations to improve NFA value. Only if the
    rectangle is not meaningful (i.e., log_nfa <= log_eps).
 */
double rect_improveORTH( struct rect * rec, image_double angles,
                            double logNT, double log_eps,
                            double* mnfa,double* mnfa_2,double* mnfa_4,
                            double*mnfap,double*mnfap_2,double*mnfap_4,
                            int Nnfa,int minsize, int minsize2,int minsize4, int orth )
{
  struct rect r;
  double log_nfa,log_nfa_new;
  double delta = 0.5;
  double delta_2 = delta / 2.0;
  int n;
  log_nfa=-1;
  //geometric width reduction
  double factor = 1./sqrt(2.);
  rect_copy(rec,&r);
  //MNFA markov table flag
  double p1check = (double)r.p+.0001;
  double p2check = (double)r.p/2. + 0.0001;



  log_nfa=rect_improve_update(r,angles,logNT,Nnfa,
		mnfa,mnfap,minsize,mnfa_2,mnfap_2,minsize2,mnfa_4,mnfap_4,minsize4,
    p1check,p2check,rec,log_nfa,orth);
  if( log_nfa > log_eps ) return log_nfa;

  /* try finer precisions */
  rect_copy(rec,&r);
  for(n=0; n<1; n++)
  {
    r.p /= 2.0;
    r.prec = r.p * M_PI;
    log_nfa=rect_improve_update(r,angles,logNT,Nnfa,
		  mnfa,mnfap,minsize,mnfa_2,mnfap_2,minsize2,mnfa_4,mnfap_4,minsize4,
      p1check,p2check,rec,log_nfa,orth);
  }
  if( log_nfa > log_eps ) return log_nfa;

  /* try to reduce width */
  rect_copy(rec,&r);
  for(n=0; n<5; n++)
  {
    if((r.width*factor) >= 0.5) //( (r.width - delta) >= 0.5 )
    {
      r.width *= factor; //-= delta;
      log_nfa=rect_improve_update(r,angles,logNT,Nnfa,
        mnfa,mnfap,minsize,mnfa_2,mnfap_2,minsize2,mnfa_4,mnfap_4,minsize4,
        p1check,p2check,rec,log_nfa,orth);
    }
  }
  if( log_nfa > log_eps ) return log_nfa;

  /* try to reduce one side of the rectangle */
  rect_copy(rec,&r);
  for(n=0; n<5; n++)
  {
    if((r.width*factor) >= 0.5) //( (r.width - delta) >= 0.5 )
    {
      delta_2 = r.width*(1.0-factor)/2.;
      r.x1 += -r.dy * delta_2;
      r.y1 +=  r.dx * delta_2;
      r.x2 += -r.dy * delta_2;
      r.y2 +=  r.dx * delta_2;
      //r.width -= delta;
      r.width *= factor;
      log_nfa=rect_improve_update(r,angles,logNT,Nnfa,
        mnfa,mnfap,minsize,mnfa_2,mnfap_2,minsize2,mnfa_4,mnfap_4,minsize4,
        p1check,p2check,rec,log_nfa,orth);
    }
  }
  if( log_nfa > log_eps ) return log_nfa;

  /* try to reduce the other side of the rectangle */
  rect_copy(rec,&r);
  for(n=0; n<5; n++)
  {
    if((r.width*factor) >= 0.5) //( (r.width - delta) >= 0.5 )
    {
      delta_2 = r.width*(1.0-factor)/2.;
      r.x1 -= -r.dy * delta_2;
      r.y1 -=  r.dx * delta_2;
      r.x2 -= -r.dy * delta_2;
      r.y2 -=  r.dx * delta_2;
      //r.width -= delta;
      r.width *= factor;
      log_nfa=rect_improve_update(r,angles,logNT,Nnfa,
        mnfa,mnfap,minsize,mnfa_2,mnfap_2,minsize2,mnfa_4,mnfap_4,minsize4,
        p1check,p2check,rec,log_nfa,orth);
    }
  }
  if( log_nfa > log_eps ) return log_nfa;

  /* try even finer precisions */
  rect_copy(rec,&r);
  for(n=0; n<1; n++)
  {
    r.p /= 2.0;
    r.prec = r.p * M_PI;
    log_nfa=rect_improve_update(r,angles,logNT,Nnfa,
        mnfa,mnfap,minsize,mnfa_2,mnfap_2,minsize2,mnfa_4,mnfap_4,minsize4,
        p1check,p2check,rec,log_nfa,orth);
  }
  return log_nfa;
}

/*----------------------------------------------------------------------------*/
/** Reduce the region size, by elimination the points far from the
    starting point, until that leads to rectangle with the right
    density of region points or to discard the region if too small.
 */
int reduce_region_radius( struct point * reg, int * reg_size,
                                 image_double modgrad, double reg_angle,
                                 double prec, double p, struct rect * rec,
                                 image_char used, image_double angles,
                                 double density_th )
{
  double density,rad1,rad2,rad,xc,yc;
  int i;

  /* check parameters */
  if( reg == NULL ) error("reduce_region_radius: invalid pointer 'reg'.");
  if( reg_size == NULL )
    error("reduce_region_radius: invalid pointer 'reg_size'.");
  if( prec < 0.0 ) error("reduce_region_radius: 'prec' must be positive.");
  if( rec == NULL ) error("reduce_region_radius: invalid pointer 'rec'.");
  if( used == NULL || used->data == NULL )
    error("reduce_region_radius: invalid image 'used'.");
  if( angles == NULL || angles->data == NULL )
    error("reduce_region_radius: invalid image 'angles'.");

  /* compute region points density */
  density = (double) *reg_size /
                         ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );

  /* if the density criterion is satisfied there is nothing to do */
  if( density >= density_th ) return TRUE;

  /* compute region's radius */
  xc = (double) reg[0].x;
  yc = (double) reg[0].y;
  rad1 = dist( xc, yc, rec->x1, rec->y1 );
  rad2 = dist( xc, yc, rec->x2, rec->y2 );
  rad = rad1 > rad2 ? rad1 : rad2;

  /* while the density criterion is not satisfied, remove farther pixels */
  while( density < density_th )
    {
      rad *= 0.75; /* reduce region's radius to 75% of its value */

      /* remove points from the region and update 'used' map */
      for(i=0; i<*reg_size; i++)
        if( dist( xc, yc, (double) reg[i].x, (double) reg[i].y ) > rad )
          {
            /* point not kept, mark it as NOTUSED */
            used->data[ reg[i].x + reg[i].y * used->xsize ] = NOTUSED;
            /* remove point from the region */
            reg[i].x = reg[*reg_size-1].x; /* if i==*reg_size-1 copy itself */
            reg[i].y = reg[*reg_size-1].y;
            --(*reg_size);
            --i; /* to avoid skipping one point */
          }

      /* reject if the region is too small.
         2 is the minimal region size for 'region2rect' to work. */
      if( *reg_size < 2 ) return FALSE;

      /* re-compute rectangle */
      region2rect(reg,*reg_size,modgrad,reg_angle,prec,p,rec);

      /* re-compute region points density */
      density = (double) *reg_size /
                         ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );
    }

  /* if this point is reached, the density criterion is satisfied */
  return TRUE;
}

int reduce_region_radiusORTH( struct point * reg, int * reg_size,
                                 image_double modgrad, double reg_angle,
                                 double prec, double p, struct rect * rec,
                                 image_char used, image_double angles,
                                 double density_th )
{
  double density,rad1,rad2,rad,xc,yc;
  int i;

  /* check parameters */
  if( reg == NULL ) error("reduce_region_radius: invalid pointer 'reg'.");
  if( reg_size == NULL )
    error("reduce_region_radius: invalid pointer 'reg_size'.");
  if( prec < 0.0 ) error("reduce_region_radius: 'prec' must be positive.");
  if( rec == NULL ) error("reduce_region_radius: invalid pointer 'rec'.");
  if( used == NULL || used->data == NULL )
    error("reduce_region_radius: invalid image 'used'.");
  if( angles == NULL || angles->data == NULL )
    error("reduce_region_radius: invalid image 'angles'.");

  /* compute region points density */
  density = (double) *reg_size /
                         ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );

  /* if the density criterion is satisfied there is nothing to do */
  if( density >= density_th ) return TRUE;

  /* compute region's radius */
  xc = (double) reg[0].x;
  yc = (double) reg[0].y;
  rad1 = dist( xc, yc, rec->x1, rec->y1 );
  rad2 = dist( xc, yc, rec->x2, rec->y2 );
  rad = rad1 > rad2 ? rad1 : rad2;

  /* while the density criterion is not satisfied, remove farther pixels */
  while( density < density_th )
    {
      rad *= 0.75; /* reduce region's radius to 75% of its value */

      /* remove points from the region and update 'used' map */
      for(i=0; i<*reg_size; i++)
        if( dist( xc, yc, (double) reg[i].x, (double) reg[i].y ) > rad )
          {
            /* point not kept, mark it as NOTUSED */
            used->data[ reg[i].x + reg[i].y * used->xsize ] = NOTUSED;
            /* remove point from the region */
            reg[i].x = reg[*reg_size-1].x; /* if i==*reg_size-1 copy itself */
            reg[i].y = reg[*reg_size-1].y;
            --(*reg_size);
            --i; /* to avoid skipping one point */
          }

      /* reject if the region is too small.
         2 is the minimal region size for 'region2rect' to work. */
      if( *reg_size < 2 ) return FALSE;

      /* re-compute rectangle */
      region2rectORTH(reg,*reg_size,modgrad,reg_angle,prec,p,rec);

      /* re-compute region points density */
      density = (double) *reg_size /
                         ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );
    }

  /* if this point is reached, the density criterion is satisfied */
  return TRUE;
}

/*----------------------------------------------------------------------------*/
/** Refine a rectangle.

    For that, an estimation of the angle tolerance is performed by the
    standard deviation of the angle at points near the region's
    starting point. Then, a new region is grown starting from the same
    point, but using the estimated angle tolerance. If this fails to
    produce a rectangle with the right density of region points,
    'reduce_region_radius' is called to try to satisfy this condition.
 */
int refine( struct point * reg, int * reg_size, image_double modgrad,
                   double reg_angle, double prec, double p, struct rect * rec,
                   image_char used, image_double angles, double density_th )
{
  double angle,ang_d,mean_angle,tau,density,xc,yc,ang_c,sum,s_sum;
  int i,n;

  /* check parameters */
  if( reg == NULL ) error("refine: invalid pointer 'reg'.");
  if( reg_size == NULL ) error("refine: invalid pointer 'reg_size'.");
  if( prec < 0.0 ) error("refine: 'prec' must be positive.");
  if( rec == NULL ) error("refine: invalid pointer 'rec'.");
  if( used == NULL || used->data == NULL )
    error("refine: invalid image 'used'.");
  if( angles == NULL || angles->data == NULL )
    error("refine: invalid image 'angles'.");

  /* compute region points density */
  density = (double) *reg_size /
                         ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );

  if (dist(rec->x1,rec->y1,rec->x2,rec->y2)<=rec->width) return FALSE;
  /* if the density criterion is satisfied there is nothing to do */
  if( density >= density_th ) return TRUE;

  /*------ First try: reduce angle tolerance ------*/

  /* compute the new mean angle and tolerance */
  xc = (double) reg[0].x;
  yc = (double) reg[0].y;
  ang_c = angles->data[ reg[0].x + reg[0].y * angles->xsize ];
  sum = s_sum = 0.0;
  n = 0;
  for(i=0; i<*reg_size; i++)
    {
      used->data[ reg[i].x + reg[i].y * used->xsize ] = NOTUSED;
      if( dist( xc, yc, (double) reg[i].x, (double) reg[i].y ) < rec->width )
        {
          angle = angles->data[ reg[i].x + reg[i].y * angles->xsize ];
          ang_d = angle_diff_signed(angle,ang_c);
          sum += ang_d;
          s_sum += ang_d * ang_d;
          ++n;
        }
    }
  mean_angle = sum / (double) n;
  tau = 2.0 * sqrt( (s_sum - 2.0 * mean_angle * sum) / (double) n
                         + mean_angle*mean_angle ); /* 2 * standard deviation */

  /* find a new region from the same starting point and new angle tolerance */
  tau=prec/2.0;
  region_grow(reg[0].x,reg[0].y,angles,reg,reg_size,&reg_angle,used,tau);
  
  /*prec = M_PI * ang_th / 180.0;
  p = ang_th / 180.0;*/
  p=tau/M_PI;
  /* if the region is too small, reject */
  if( *reg_size < 2 ) return FALSE;

  /* re-compute rectangle */
  region2rect(reg,*reg_size,modgrad,reg_angle,tau,p,rec);

  /* re-compute region points density */
  density = (double) *reg_size /
                      ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );

  /*------ Second try: reduce region radius ------*/
  if( density < density_th )
    return reduce_region_radius( reg, reg_size, modgrad, reg_angle, prec, p,
                                 rec, used, angles, density_th );


  if (dist(rec->x1,rec->y1,rec->x2,rec->y2)<=rec->width) return FALSE;
  /* if this point is reached, the density criterion is satisfied */
  return TRUE;
}

/*----------------------------------------------------------------------------*/
/** Refine a rectangle.

    For that, an estimation of the angle tolerance is performed by the
    standard deviation of the angle at points near the region's
    starting point. Then, a new region is grown starting from the same
    point, but using the estimated angle tolerance. If this fails to
    produce a rectangle with the right density of region points,
    'reduce_region_radius' is called to try to satisfy this condition.
 */
int refineORTH( struct point * reg, int * reg_size, image_double modgrad,
                   double reg_angle, double prec, double p, struct rect * rec,
                   image_char used, image_double angles, double density_th ,int orth)
{
  double angle,ang_d,mean_angle,tau,density,xc,yc,ang_c,sum,s_sum;
  int i,n;

  /* check parameters */
  if( reg == NULL ) error("refine: invalid pointer 'reg'.");
  if( reg_size == NULL ) error("refine: invalid pointer 'reg_size'.");
  if( prec < 0.0 ) error("refine: 'prec' must be positive.");
  if( rec == NULL ) error("refine: invalid pointer 'rec'.");
  if( used == NULL || used->data == NULL )
    error("refine: invalid image 'used'.");
  if( angles == NULL || angles->data == NULL )
    error("refine: invalid image 'angles'.");
  if((orth!=1)&&(orth!=0)) error("refineORTH: orth flag not 0/1");

  /* compute region points density */
  density = (double) *reg_size /
                         ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );

  /* if the density criterion is satisfied there is nothing to do */
  if (dist(rec->x1,rec->y1,rec->x2,rec->y2)<=rec->width) return FALSE;
  if( density >= density_th ) return TRUE;

  /*------ First try: reduce angle tolerance ------*/

  /* compute the new mean angle and tolerance */
  xc = (double) reg[0].x;
  yc = (double) reg[0].y;
  ang_c = angles->data[ reg[0].x + reg[0].y * angles->xsize ];
  sum = s_sum = 0.0;
  n = 0;
  for(i=0; i<*reg_size; i++)
    {
      used->data[ reg[i].x + reg[i].y * used->xsize ] = NOTUSED;
      if( dist( xc, yc, (double) reg[i].x, (double) reg[i].y ) < rec->width )
        {
          angle = angles->data[ reg[i].x + reg[i].y * angles->xsize ];
          ang_d = angle_diff_signed(angle,ang_c);
          sum += ang_d;
          s_sum += ang_d * ang_d;
          ++n;
        }
    }
  mean_angle = sum / (double) n;
  tau = 2.0 * sqrt( (s_sum - 2.0 * mean_angle * sum) / (double) n
                         + mean_angle*mean_angle ); /* 2 * standard deviation */

  /* find a new region from the same starting point and new angle tolerance */
  tau=prec/2.0;
  region_growORTH( reg[0].x,reg[0].y, angles, 
                      reg,reg_size, &reg_angle, &reg_angle, used, tau);

  /*prec = M_PI * ang_th / 180.0;
  p = ang_th / 180.0;*/
  p=tau/M_PI;
  /* if the region is too small, reject */
  if( *reg_size < 2 ) return FALSE;

  /* re-compute rectangle */
  region2rectORTH(reg,*reg_size,modgrad,reg_angle,tau,p,rec);
 


  /* re-compute region points density */
  density = (double) *reg_size /
                      ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );

  /*------ Second try: reduce region radius ------*/
  if( density < density_th )
    return reduce_region_radiusORTH( reg, reg_size, modgrad, reg_angle, prec, p,
                                 rec, used, angles, density_th );

  
  if (dist(rec->x1,rec->y1,rec->x2,rec->y2)<=rec->width) return FALSE;
  /* if this point is reached, the density criterion is satisfied */
  return TRUE;
}
