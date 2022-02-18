/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/gradient.c: compute gradient image using Sobel or Gradient-by-Ratio in 2D/3D

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
/*---------------------------------- Import ---------------------------------*/
/*----------------------------------------------------------------------------*/

#ifndef GSOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include <tgmath.h>
#include <limits.h>
#include <float.h>
#include<string.h>
#include <time.h>
#include<gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>
#include <sys/mman.h>

#include "gradient.h"
#include "constants.h"
#include "misc.h"
#include "tuples.h"

/*----------------------------------------------------------------------------*/
/*------------------------------ 2D Gradient ---------------------------------*/
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Computes the direction of the level line of 'in' at each point.

    The result is:
    - an image_double with the angle at each pixel, or NOTDEF if not defined.
    - the image_double 'modgrad' (a pointer is passed as argument)
      with the gradient magnitude at each point.
    - a list of pixels 'list_p' roughly ordered by decreasing
      gradient magnitude. (The order is made by classifying points
      into bins by gradient magnitude. The parameters 'n_bins' and
      'max_grad' specify the number of bins and the gradient modulus
      at the highest bin. The pixels in the list would be in
      decreasing gradient magnitude, up to a precision of the size of
      the bins.)
    - a pointer 'mem_p' to the memory used by 'list_p' to be able to
      free the memory when it is not used anymore.
 */
image_double ll_angle( image_double in,
                              struct coorlist ** list_p, void ** mem_p,
                              image_double * modgrad, unsigned int n_bins,double alpha)
{
  image_double g;
  double norm;
  unsigned int n,p,x,y,adr,i; 
  int list_count = 0;
  struct coorlist * list;
  struct coorlist ** range_l_s; /* array of pointers to start of bin list */
  struct coorlist ** range_l_e; /* array of pointers to end of bin list */
  struct coorlist * start;
  struct coorlist * end;
  double max_grad = 0.0;

  /* check parameters */
  if( in == NULL || in->data == NULL || in->xsize == 0 || in->ysize == 0 )
    error("ll_angle: invalid image.");
 
  if( list_p == NULL ) error("ll_angle: NULL pointer 'list_p'.");
  if( mem_p == NULL ) error("ll_angle: NULL pointer 'mem_p'.");
  if( modgrad == NULL ) error("ll_angle: NULL pointer 'modgrad'.");
  if( n_bins == 0 ) error("ll_angle: 'n_bins' must be positive.");

  /* image size shortcuts */
  n = in->ysize;
  p = in->xsize;
  /* allocate output image */
  g = new_image_double(in->xsize,in->ysize);

  /* get memory for the image of gradient modulus */
  *modgrad = new_image_double(in->xsize,in->ysize);

  /* get memory for "ordered" list of pixels */
  list = (struct coorlist *) calloc( (size_t) (n*p), sizeof(struct coorlist) );
  *mem_p = (void *) list;
  range_l_s = (struct coorlist **) calloc( (size_t) n_bins,
                                           sizeof(struct coorlist *) );
  range_l_e = (struct coorlist **) calloc( (size_t) n_bins,
                                           sizeof(struct coorlist *) );
  if( list == NULL || range_l_s == NULL || range_l_e == NULL )
    error("not enough memory.");
  for(i=0;i<n_bins;i++) range_l_s[i] = range_l_e[i] = NULL;



  /* 'undefined' on the down and righ t boundaries */
  for(x=0;x<p;x++) g->data[(n-1)*p+x] = NOTDEF;
  for(y=0;y<n;y++) g->data[p*y+p-1]   = NOTDEF;

  
  int j,k;
  double ax,ay,an,ap;
  int X,Y;
  Y=n;  X=p;
  double grfilt = 0;
  if (alpha<0)
  {
    grfilt = -alpha;
    alpha = 0.;
  }
  double beta=alpha;
  int nb=16;

  /* create a simple image: left half black, right half gray */
  double * gradx;
  gradx= (double *) malloc( X * Y * sizeof(double) );
  double * grady;
  grady= (double *) malloc( X * Y * sizeof(double) );
  int imgSize=X*Y;
  int bitdepth=16;
  double *img1;
  img1= (double *) malloc( X * Y * sizeof(double) );
  double *img2;
  img2= (double *) malloc( X * Y * sizeof(double) );
  for(i=0;i<imgSize;i++)
  {
    /*img1[i]=pow(image[i],2);*/
    img1[i]=in->data[i];
    if(img1[i]<1.) img2[i]=1.;
    else img2[i]=img1[i];
  }
  int longueur=ceil(log(10)*beta);
  /*longueur=wid;*/
  int sobel=0;
  double sobel1=5./8.;//10.;//2
  double sobel2=3./16.;//3.;//1
  if (longueur==0)
  {
    longueur=1;
    sobel=1;
    
  }
  int largeur=longueur;
  double * gradx1;
  gradx1= (double *) malloc( X * Y * sizeof(double) );
  double * grady1;
  grady1= (double *) malloc( X * Y * sizeof(double) );
  for(j=0;j<Y;j++)
  {
    for(i=0;i<X;i++)
    {
      double Mx=0.;
      double My=0.;
      for(k=-largeur;k<=largeur;k++)
      {
        int xk=min1(max1(i+k,0),X-1);
        int yk=min1(max1(j+k,0),Y-1);
        if(sobel>=1)
        {
          Mx += k*img1[xk+j*X];  
          if (k==0) {
            My += sobel1*img1[xk+j*X];  
          }
          else {
            My += sobel2*img1[xk+j*X];  
          }                         
        }
        else
        {
          double coeff=exp(-(double) abs(k)/beta);
          Mx+=coeff*img2[xk+j*X];
          My+=coeff*img2[i+yk*X];
        }
      }
      gradx1[i+j*X]=Mx;
      grady1[i+j*X]=My;
    }
  }
  for(j=0;j<Y;j++)
  {
    for(i=0;i<X;i++)
    {
      double Mxg=0;
      double Mxd=0;
      double Myg=0;
      double Myd=0;
      if(sobel>=1)
      {
        int yl1;
        if(j-1<0) yl1=0;
        else yl1=j-1;
        int yl2;
        if(j+1>Y-1)  yl2=Y-1;
        else yl2=j+1;
        
        Myg += grady1[i+yl2*X];  
        Myg -= grady1[i+yl1*X];  
        Mxg += sobel1*gradx1[i+j*X];  
        Mxg += sobel2*gradx1[i+yl2*X];  
        Mxg += sobel2*gradx1[i+yl1*X];  
      }
      else
      {
        for(k=1;k<=longueur;k++)
        {
          double coeff=exp(-(double) abs(k)/beta);
          int yl1;
          if(j-k<0) yl1=0;
          else yl1=j-k;
          int yl2;
          if(j+k>Y-1)  yl2=Y-1;
          else yl2=j+k;
          
          Mxg+=coeff*gradx1[i+yl1*X];
          Mxd+=coeff*gradx1[i+yl2*X];
          int xl1=max1(i-k,0);
          int xl2=min1(i+k,X-1);;
          Myg+=coeff*grady1[xl1+j*X];
          Myd+=coeff*grady1[xl2+j*X];
        }
      }
      if(sobel>=1)
      {
        gradx[i+j*X]=Mxg;
        grady[i+j*X]=Myg;
      }
      else
      {
        gradx[i+j*X]=log(Mxd/Mxg);
        grady[i+j*X]=log(Myd/Myg);
      }
    }
  }
  for(i=0;i<X;i++)
  {
    for(j=0;j<Y;j++)
    {
      adr = j*X+i;
      ay=gradx[adr];
      ax=grady[adr];
      an=(double) hypot((double) ax,(double) ay);
      norm=an;
      if (norm<=grfilt)
      {
          norm=0.;   
      }
      (*modgrad)->data[adr] = norm; /* store gradient norm */
      if( norm <= 0.0 ) /* norm too small, gradient no defined */
        g->data[adr] = NOTDEF; /* gradient angle not defined */
      else
      {
        /* gradient angle computation */
        ap=atan2((double) ax,-(double) ay);
        g->data[adr] = ap;
        /* look for the maximum of the gradient */
        if( norm > max_grad ) max_grad = norm;
      }       
    }
  }
  int i0;
  /* compute histogram of gradient values */
  for(x=0;x<X-1;x++)
    for(y=0;y<Y-1;y++)
    {
      norm = (*modgrad)->data[y*p+x];

      /* store the point in the right bin according to its norm */
      i0= (unsigned int) (norm * (double) n_bins / max_grad);
      if( i0 >= n_bins ) i0 = n_bins-1;
      if( range_l_e[i0] == NULL )
        range_l_s[i0] = range_l_e[i0] = list+list_count++;
      else
        {
          range_l_e[i0]->next = list+list_count;
          range_l_e[i0] = list+list_count++;
        }
      range_l_e[i0]->x = (int) x;
      range_l_e[i0]->y = (int) y;
      range_l_e[i0]->next = NULL;
    }

  /* Make the list of pixels (almost) ordered by norm value.
     It starts by the larger bin, so the list starts by the
     pixels with the highest gradient value. Pixels would be ordered
     by norm value, up to a precision given by max_grad/n_bins.
   */
  for(i=n_bins-1; i>0 && range_l_s[i]==NULL; i--);
  start = range_l_s[i];
  end = range_l_e[i];
  if( start != NULL )
    while(i>0)
      {
        --i;
        if( range_l_s[i] != NULL )
          {
            end->next = range_l_s[i];
            end = range_l_e[i];
          }
      }
  *list_p = start;

  /* free memory */
  free( (void *) range_l_s );
  free( (void *) range_l_e );
  free ((void *) gradx);
  free((void *) grady);
  free((void *) img1);
  free((void *) img2);
  free((void *) gradx1);
  free((void *) grady1);
  return g;
}


/*----------------------------------------------------------------------------*/
/*------------------------------ 3D Gradient ---------------------------------*/
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/** Create a new image3_double pair to store gradient orientations. */
grads new_grads(unsigned int xsize, unsigned int ysize, unsigned int zsize)
{
  grads image;
  /* check parameters */
  if(xsize==0 || ysize==0 || zsize==0) error("new_image3_double: invalid image size.");
  /* get memory */
  image = (grads) malloc( sizeof(struct grads_s) );
  if( image == NULL ) error("not enough memory.");

  image->az = new_image3_double(xsize,ysize,zsize);
  image->el = new_image3_double(xsize,ysize,zsize); 
  if( image->az == NULL || image->az->data == NULL ) error("not enough memory.");
  if( image->el == NULL || image->el->data == NULL ) error("not enough memory.");
  return image;
}

/*----------------------------------------------------------------------------*/
/** Free memory used in gradient structure. */
void free_grads(grads i)
{
  if( i == NULL || i->az->data == NULL )
    error("free_grads: invalid input image.");
  free_image3_double(i->az );
  free_image3_double(i->el );
  free( (void *) i );
}

/*----------------------------------------------------------------------------*/
/** Computes the direction of the GRADIENT VECTOR 
 * of 'in' at each point, in constast to the 2D alg., 
 * since the 'Level Line' is an ambiguous 2D surface in 3D.
 * The GR algorithm is identical to 2D, but the atan is oriented normally.  
 * This correction mandates the following changes:
 *
 * get_theta3:   Comuputes minimum eigenvector for principal axis in region2rect3
 * region3_grow: Grow regions with parallel neighbors 
 * rect3_nfa:    WOULD LIKE to check orthogonality alignment to principle axis, 
 *                however to preserve parallel statistics, we instead count 
 *                parallel alignments to each principal axis and count the maximum.
 *                TODO: test switch to orthogonality counting
 * make_markov3: Count paralle contributions to horiz/vert/depth lines
 *
 * The alg. exploits  seperability for fast runtime, 
 * and is accelerated by MPI parallelization on the outermost loop. 
 *
 *   The result is:
 *   - an image_double with the azimuth at each pixel, or NOTDEF if not defined.
 *   - an image_double with the elevation at each pixel, or NOTDEF if not defined.
 *   - the image_double 'modgrad' (a pointer is passed as argument)
 *     with the gradient magnitude at each point.
 *   - a list of pixels 'list_p' roughly ordered by decreasing
 *     gradient magnitude. (The order is made by classifying points
 *     into bins by gradient magnitude. The parameters 'n_bins' and
 *     'max_grad' specify the number of bins and the gradient modulus
 *     at the highest bin. The pixels in the list would be in
 *     decreasing gradient magnitude, up to a precision of the size of
 *     the bins.)
 *   - a pointer 'mem_p' to the memory used by 'list_p' to be able to
 *     free the memory when it is not used anymore.
 */
grads ll_angle3( image3_double in,
                        struct coorlist3 ** list_p, void ** mem_p,
                        image3_double * modgrad, 
                        unsigned int n_bins,double alpha)
{
  unsigned int m,n,p,x,y,z,adr,i;
  double norm;
  int list_count = 0;
  struct coorlist3 * list;
  struct coorlist3 ** range_l_s; /* array of pointers to start of bin list */
  struct coorlist3 ** range_l_e; /* array of pointers to end of bin list */
  struct coorlist3 * start;
  struct coorlist3 * end;
  double max_grad = 0.0;

  /* check parameters */
  if( in == NULL || in->data == NULL || in->xsize == 0 || in->ysize == 0 || in->zsize ==0 )
    error("ll_angle3: invalid image."); 
  if( list_p == NULL ) error("ll_angle3: NULL pointer 'list_p'.");
  if( mem_p == NULL ) error("ll_angle3: NULL pointer 'mem_p'.");
  if( modgrad == NULL ) error("ll_angle3: NULL pointer 'modgrad'.");
  if( n_bins == 0 ) error("ll_angle3: 'n_bins' must be positive.");

  /* image size shortcuts */
  m = in->zsize;
  n = in->ysize;
  p = in->xsize;

  /* get memory for the gradient angle data */
  grads angles = new_grads(p,n,m);
  /* get memory for the image of gradient modulus */
  *modgrad = new_image3_double(p,n,m);
  /* get memory for "ordered" list of pixels */
  list = (struct coorlist3 *) calloc( (size_t) (m*n*p), sizeof(struct coorlist3) );
  *mem_p = (void *) list;
  range_l_s = (struct coorlist3 **) calloc( (size_t) n_bins,
                                           sizeof(struct coorlist *) );
  range_l_e = (struct coorlist3 **) calloc( (size_t) n_bins,
                                           sizeof(struct coorlist *) );
  if( list == NULL || range_l_s == NULL || range_l_e == NULL )
    error("ll_angle3: not enough memory in list or range.");
  for(i=0;i<n_bins;i++) range_l_s[i] = range_l_e[i] = NULL;
  
  //undefined outermost boundaries
  /* 'undefined' on the down and right (outermost)  boundaries 
   * undefined on the z-axis/depth boundaries, aka 'outer'
   */
  //Let y=(n-1) for 'down boundary'.  Cycle all x and z
  for(x=0;x<p;x++) for(z=0;z<m;z++){angles->az->data[ (z) + m*( (x) + p*(n-1) ) ] = NOTDEF;
                    angles->el->data[ (z) + m*( (x) + p*(n-1) ) ] = NOTDEF;}
  //Let x=(p-1), for 'right boundary'.  Cycle all y and z
  for(y=0;y<n;y++) for(z=0;z<m;z++){angles->az->data[ (z) + m*( (p-1) + p*(y) ) ] = NOTDEF;
                    angles->el->data[ (z) + m*( (p-1) + p*(y) ) ] = NOTDEF;}
  //Let z=(m-1), for 'out boundary'.  Cycle all x and z
  for(x=0;x<p;x++) for(y=0;y<n;y++){angles->az->data[ (m-1) + m*( (x) + p*(y) ) ] = NOTDEF;
                    angles->el->data[ (m-1) + m*( (x) + p*(y) ) ] = NOTDEF;}
 

  int j,k,ha,hb;
  double ax,ay,az,an;
  int X,Y,Z;
  Z=m;Y=n;X=p;
  int imgSize=X*Y*Z;
  double grfilt = 0;
  if (alpha<0)
  {
    grfilt = -alpha;
    alpha = 0.;
  }
  double beta=alpha;
  int nb=16;
  int bitdepth=16;
  double imgtemp;
  int longueur=ceil(log(10)*beta);
  int sobel=0;
  double sobel1=5./8.;//10.;//2
  double sobel2=3./16.;//3.;//1
  if (longueur==0)
  {
    longueur=1;
    sobel=1;
  }
  int largeur=longueur;
  double Mg, Md, coeff;
  double Mgx,Mdx,Mgy,Mdy,Mgz,Mdz;
  int h,xx,yy,zz,hx,hy,hz;

  //printf("Get Gradients\n");fflush(stdout);
  /*Temporary memory matrices for re-use of kernel data*/
  double *img1=   (double *) malloc( imgSize *  sizeof(double) );
  if (img1==NULL) error("not enough memory.");
  double *gradx1= (double *) malloc( imgSize *  sizeof(double) );
  if (gradx1==NULL) error("not enough memory.");
  double *gradx2= (double *) malloc( imgSize *  sizeof(double) );
  if (gradx2==NULL) error("not enough memory.");
  double *grady1= (double *) malloc( imgSize *  sizeof(double) );
  if (grady1==NULL) error("not enough memory.");
  double *grady2= (double *) malloc( imgSize *  sizeof(double) );
  if (grady2==NULL) error("not enough memory.");
  double *gradz1= (double *) malloc( imgSize *  sizeof(double) );
  if (gradz1==NULL) error("not enough memory.");
  double *gradz2= (double *) malloc( imgSize *  sizeof(double) );
  if (gradz2==NULL) error("not enough memory.");


  for(i=0;i<imgSize;i++)
  {
    //required data clipping
    img1[i]=(double)in->data[i];
    if(img1[i]<1.) img1[i]=1.;
    //fill temp matrices to ensure space 
    gradx1[i]=(double)in->data[i];
    gradx2[i]=(double)in->data[i];
    grady1[i]=(double)in->data[i];
    grady2[i]=(double)in->data[i];
    gradz1[i]=(double)in->data[i];
    gradz2[i]=(double)in->data[i];
  }

  /*Wall-clock timing */
  double startT,endT;
  //int zcheck = largeur-1;
  double zbeta = beta;//-0.5;//.5;//.5;//.25;
  int xyflag = 0;
  int zcheck=ceil(log(10)*zbeta);
  startT=omp_get_wtime();
  #pragma omp parallel default(none) shared(sobel1,sobel2,xyflag,zbeta,zcheck,sobel,X,Y,Z,largeur,beta,img1,gradx1,grady1,gradz1) private(i,j,k,h,xx,yy,zz,coeff,Mdx,Mdy,Mdz)   
  {
    #pragma omp for 
    for(k=0;k<Z;k++){
      for(j=0;j<Y;j++) {
        for(i=0;i<X;i++) {
          Mdx=0.;  Mdy=0.;  Mdz=0.; 
          for(h=-largeur;h<=largeur;h++)
          {     
            xx = (i+h)>0?(i+h):0; xx = xx<(X-1)?xx:(X-1);
            yy = (j+h)>0?(j+h):0; yy = yy<(Y-1)?yy:(Y-1);
            zz = (k+h)>0?(k+h):0; zz = zz<(Z-1)?zz:(Z-1);
            if(sobel>=1)
            {
              Mdx += h*img1[k+Z*(i+X*yy)];  
              if (h==0) {
                Mdy += sobel1*img1[k+Z*(i+X*yy)];  
                Mdz += sobel1*img1[k+Z*(i+X*yy)];  
                //Mdy += 2.*img1[zz+Z*(i+X*j)];  
                //Mdz += 2.*img1[k+Z*(xx+X*j)];  
              }
              else {
                Mdy += sobel2*img1[k+Z*(i+X*yy)];  
                Mdz += sobel2*img1[k+Z*(i+X*yy)];  
                //Mdy += 1.*img1[zz+Z*(i+X*j)];  
                //Mdz += 1.*img1[k+Z*(xx+X*j)];  
              }                         
            }
            else
            {
              coeff=exp(-(double) abs(h)/beta);
              if (Z>1)
              {  
                if (xyflag==0)
                {
                  Mdx+=coeff*img1[k +Z*(i +X*yy)];
                  if(abs(h)<=zcheck) //for 5x5 kernel, max at 3
                  {Mdy+=exp(-(double)abs(h)/zbeta)*img1[zz+Z*(i +X*j )];}
                  Mdz+=coeff*img1[k +Z*(xx+X*j )];
                }
                else
                {
                  Mdy+=coeff*img1[zz+Z*(i +X*j )];
                  if(abs(h)<=zcheck) //for 5x5 kernel, max at 3
                  {
                    Mdx+=exp(-(double)abs(h)/zbeta)*img1[k +Z*(i +X*yy)];
                    Mdz+=exp(-(double)abs(h)/zbeta)*img1[k +Z*(xx+X*j )];
                  }
                }

              }
              else
              {
                Mdx+=coeff*img1[k +Z*(i +X*yy)];
                Mdy = img1[k+Z*(i+X*j)];
              }
            }
          }
          gradx1[k+Z*(i+X*j)]=Mdx; grady1[k+Z*(i+X*j)]=Mdy; gradz1[k+Z*(i+X*j)]=Mdz;
  }}}}
  endT=omp_get_wtime();
  startT = omp_get_wtime();
  #pragma omp parallel default(none) shared(sobel1,sobel2,xyflag,zbeta,sobel,zcheck,X,Y,Z,largeur,beta,gradx1,grady1,gradz1,gradx2,grady2,gradz2) private(i,j,k,h,xx,yy,zz,coeff,Mdx,Mdy,Mdz)   
  {
    #pragma omp for 
    for(k=0;k<Z;k++){
      for(j=0;j<Y;j++) {
        for(i=0;i<X;i++) {
          Mdx=0.;  Mdy=0.;  Mdz=0.; 
          for(h=-largeur;h<=largeur;h++)
          {     
            xx = (i+h)>0?(i+h):0; xx = xx<(X-1)?xx:(X-1);
            yy = (j+h)>0?(j+h):0; yy = yy<(Y-1)?yy:(Y-1);
            zz = (k+h)>0?(k+h):0; zz = zz<(Z-1)?zz:(Z-1);
            if(sobel>=1)
            {
              Mdy += h*grady1[k+Z*(xx+X*j)];  
              if (h==0)
              {
                Mdx += sobel1*gradx1[k+Z*(xx+X*j)];  
                Mdz += sobel1*gradz1[k+Z*(xx+X*j)]; 
                //Mdx += 2.*gradx1[zz+Z*(i+X*j)];  
                //Mdz += 2.*gradz1[k+Z*(i+X*yy)]; 
              }
              else
              {
                Mdx += sobel2*gradx1[k+Z*(xx+X*j)];  
                Mdz += sobel2*gradz1[k+Z*(xx+X*j)]; 
                //Mdx += 1.*gradx1[zz+Z*(i+X*j)];  
                //Mdz += 1.*gradz1[k+Z*(i+X*yy)]; 
              }                         
            }
            else
            {
              coeff=exp(-(double) abs(h)/beta);
              if (Z>1)
              {  
                if (xyflag==0)
                {
                  if(abs(h)<=zcheck) //for 5x5 kernel, max at 3
                  {Mdx+=exp(-(double)abs(h)/zbeta)*gradx1[zz+Z*(i +X*j)];}
                  Mdy+=coeff*grady1[k +Z*(xx+X*j )];
                  Mdz+=coeff*gradz1[k +Z*(i +X*yy)];
                }
                else
                {
                  Mdx+=coeff*gradx1[zz+Z*(i +X*j)];
                  if(abs(h)<=zcheck) //for 5x5 kernel, max at 3
                  {
                    Mdy+=exp(-(double)abs(h)/zbeta)*grady1[k +Z*(xx +X*j)];
                    Mdz+=exp(-(double)abs(h)/zbeta)*gradz1[k +Z*(i+X*yy )];
                  }
                } 
              }
              else
              {
                Mdx = gradx1[k+Z*(i+X*j)];
                Mdy+=coeff*grady1[k +Z*(xx+X*j )];
              }
            }
          }
          gradx2[k+Z*(i+X*j)]=Mdx;
          grady2[k+Z*(i+X*j)]=Mdy;
          gradz2[k+Z*(i+X*j)]=Mdz;
  }}}}
  endT=omp_get_wtime();

  //xyflag=(int)(xyflag==0);
  double an2;
  #pragma omp parallel default(none) shared(grfilt,sobel1,sobel2,xyflag,zbeta,zcheck,sobel,X,Y,Z,largeur,beta,gradx2,grady2,gradz2,angles,modgrad) private(i,j,k,h,xx,yy,zz,coeff,Mdx,Mdy,Mdz,Mgx,Mgy,Mgz,ay,ax,az,adr,an,an2)   
  {
    #pragma omp for 
    for(k=0;k<Z;k++){
      for(j=0;j<Y;j++) {
        for(i=0;i<X;i++) {
          Mdx=0.; Mgx=0.; Mdy=0.; Mgy=0.; Mdz=0.; Mgz=0.;
          for(h=-largeur;h<=largeur;h++)
          {     
            xx = (i+h)>0?(i+h):0; xx = xx<(X-1)?xx:(X-1);
            yy = (j+h)>0?(j+h):0; yy = yy<(Y-1)?yy:(Y-1);
            zz = (k+h)>0?(k+h):0; zz = zz<(Z-1)?zz:(Z-1);
            
  
            if(sobel>=1)
            {
              Mdz += h*gradz2[zz+Z*(i+X*j)];  
              if (h==0)
              {
                Mdx += sobel1*gradx2[zz+Z*(i+X*j)];  
                Mdy += sobel1*grady2[zz+Z*(i+X*j)];  
                //Mdx += 2.*gradx2[k+Z*(xx+X*j)];  
                //Mdy += 2.*grady2[k+Z*(i+X*yy)];  
              }
              else
              {
                Mdx += sobel2*gradx2[zz+Z*(i+X*j)];  
                Mdy += sobel2*grady2[zz+Z*(i+X*j)];  
                //Mdx += 1.*gradx2[k+Z*(xx+X*j)];  
                //Mdy += 1.*grady2[k+Z*(i+X*yy)];  
              }                         
            }
            else
            {
              coeff=exp(-(double) abs(h)/beta);
              if (xyflag==0)
              {
                if(h<0)
                {
                  Mgx+=coeff*gradx2[k +Z*(xx+X*j)];
                  Mgy+=coeff*grady2[k +Z*(i +X*yy )];
                  if (Z>1) {
                    if(abs(h)<=zcheck) //for 5x5 kernel, max at 3
                    Mgz+=exp(-(double)abs(h)/zbeta)*gradz2[zz+Z*(i +X*j )];}
                }
                if(h>0)
                {
                  Mdx+=coeff*gradx2[k +Z*(xx+X*j)];
                  Mdy+=coeff*grady2[k +Z*(i +X*yy )];
                  if (Z>1) {
                    if(abs(h)<=zcheck) //for 5x5 kernel, max at 3
                    Mdz+=exp(-(double)abs(h)/zbeta)*gradz2[zz+Z*(i +X*j )];}
                }


              }
              else
              {
                if(h<0)
                {
                  Mgz+=coeff*gradz2[zz+Z*(i +X*j )];
                  if (Z>1) {
                    if(abs(h)<=zcheck) //for 5x5 kernel, max at 3
                    {
                    Mgx+=exp(-(double)abs(h)/zbeta)*gradx2[k +Z*(xx+X*j)];
                    Mgy+=exp(-(double)abs(h)/zbeta)*grady2[k +Z*(i +X*yy )];
                    }}
                }
                if(h>0)
                {
                  Mdz+=coeff*gradz2[zz+Z*(i +X*j )];
                  if (Z>1) {
                    if(abs(h)<=zcheck) //for 5x5 kernel, max at 3
                    {
                    Mdx+=exp(-(double)abs(h)/zbeta)*gradx2[k +Z*(xx+X*j)];
                    Mdy+=exp(-(double)abs(h)/zbeta)*grady2[k +Z*(i +X*yy )];
                    }}
                }
              } 
              
            }
          }
          
          //SWAP AX AND AY per LSD ordering 
          if(sobel>=1)
          {
            ay=(double)Mdx;
            ax=(double)Mdy;
            az=(double)Mdz;
          }
          else
          { 
            ay=(double)log(Mdx/Mgx);
            ax=(double)log(Mdy/Mgy);
            if (Z>1) az=(double)log(Mdz/Mgz);
            else az=0;
          }
          adr = (unsigned int)  k+Z*(j*X+i);
          an= (double)sqrt(ax*ax + ay*ay + az*az); 
          an2= (double)sqrt(ax*ax + ay*ay); 
  
          if (an<=grfilt)
          {
              an=0.;   
          }
         //Mgx=0 at (i,j)=(0,0).  
         //This is expected (see 'undefined outmost boundaries' above)
         //The 0-an setting enforces the NOTDEF angle condition
         if(!isfinite((double)an)) an = 0.0; //correction attempt
         (*modgrad)->data[adr] =  an; /* store gradient norm */

          if( an <= 0.0 ) /* norm too small, gradient no defined */
          {
            angles->az->data[adr] = NOTDEF;
            angles->el->data[adr] = NOTDEF;
          }
          else
          {
            /* gradient angle computation */
            if(an2<=0.0) {angles->az->data[adr] =0.0;}
            else {angles->az->data[adr] = atan2(ay,ax);}
            angles->el->data[adr] = acos(az/an);
          } 
  }}}}

  /*Compute maximum gradient after computation to avoid OMP race condition*/ 
  for(k=0;k<Z;k++){
    for(j=0;j<Y;j++) {
      for(i=0;i<X;i++) {
        {an= (double) (*modgrad)->data[k+Z*(j*X+i)]; 
        }//printf("%2.f",an);fflush(stdout);}
      if( an > max_grad ) max_grad = an;
  }}}
  //printf("maxgrad: %.2f\n",max_grad);fflush(stdout);

  /*free temporary data matrices*/
  /*
  free((void *) img1);
  free((void *) gradx1);
  free((void *) grady1);
  free((void *) gradz1);
  free((void *) gradx2);
  free((void *) grady2);
  free((void *) gradz2);
  */
  free(img1);
  free(gradx1);
  free(grady1);
  free(gradz1);
  free(gradx2);
  free(grady2);
  free(gradz2);

  int i0;
  /* compute histogram of gradient values */
  for(x=0;x<X-1;x++){
    for(y=0;y<Y-1;y++){
      for(z=0;z<Z-1;z++){
        norm = (*modgrad)->data[z+m*(x+y*p)];
        /* store the point in the right bin according to its norm */
        i0= (unsigned int) (norm * (double) n_bins / max_grad);
        if( i0 >= n_bins ) i0 = n_bins-1;
        if( range_l_e[i0] == NULL )
          range_l_s[i0] = range_l_e[i0] = list+list_count++;
        else
        {
          range_l_e[i0]->next = list+list_count;
          range_l_e[i0] = list+list_count++;
        }
        range_l_e[i0]->x = (int) x;
        range_l_e[i0]->y = (int) y;
        range_l_e[i0]->z = (int) z;
        range_l_e[i0]->next = NULL;
  }}}
  
  /* Make the list of pixels (almost) ordered by norm value.
     It starts by the larger bin, so the list starts by the
     pixels with the highest gradient value. Pixels would be ordered
     by norm value, up to a precision given by max_grad/n_bins.
   */
  
  for(i=n_bins-1; i>0 && range_l_s[i]==NULL; i--);
  start = range_l_s[i];
  end = range_l_e[i];
  if( start != NULL )
    while(i>0)
    {
      --i;
      if( range_l_s[i] != NULL )
      {
        end->next = range_l_s[i];
        end = range_l_e[i];
      }
    }
  *list_p = start;
  free( (void *) range_l_s );
  free( (void *) range_l_e );
  int err=0;

  //printf("LL_Grad Completed\n");fflush(stdout);
  return angles;
}

/*----------------------------------------------------------------------------*/
/*------------------------------ Alignments ---------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Is point (x,y) aligned to angle theta, up to precision 'prec'?
 */
int isaligned( int x, int y, image_double angles, double theta,
                      double prec )
{
  double a;

  /* check parameters */
  if( angles == NULL || angles->data == NULL )
    error("isaligned: invalid image 'angles'.");
  if( x < 0 || y < 0 || x >= (int) angles->xsize || y >= (int) angles->ysize )
    error("isaligned: (x,y) out of the image.");
  if( prec < 0.0 ) error("isaligned: 'prec' must be positive.");

  /* angle at pixel (x,y) */
  a = angles->data[ x + y * angles->xsize ];

  /* pixels whose level-line angle is not defined
     are considered as NON-aligned */
  if( a == NOTDEF ) return FALSE;  /* there is no need to call the function
                                      'double_equal' here because there is
                                      no risk of problems related to the
                                      comparison doubles, we are only
                                      interested in the exact NOTDEF value */

  /* it is assumed that 'theta' and 'a' are in the range [-pi,pi] */
  theta -= a;
  if( theta < 0.0 ) theta = -theta;
  if( theta > M_3_2_PI )
    {
      theta -= M_2__PI;
      if( theta < 0.0 ) theta = -theta;
    }

  return theta <= prec;
}

/*----------------------------------------------------------------------------*/
/** Is point (x,y) aligned to angle theta, up to precision 'prec'?
 */
int isalignedORTH( int x, int y, image_double angles, double theta,
                      double prec )
{
  double a;

  /* check parameters */
  if( angles == NULL || angles->data == NULL )
    error("isaligned: invalid image 'angles'.");
  if( x < 0 || y < 0 || x >= (int) angles->xsize || y >= (int) angles->ysize )
    error("isaligned: (x,y) out of the image.");
  if( prec < 0.0 ) error("isaligned: 'prec' must be positive.");
  
  double par1 = theta+0;
  if (par1>M_PI) par1-=M_PI;
  double par2 = theta-M_PI;
  if (par2<(-1.*M_PI)) par2+=M_PI;
  return (isaligned(x,y,angles,par1,prec)||isaligned(x,y,angles,par2,prec));
}






/*----------------------------------------------------------------------------*/
/** Is point (x,y,z) aligned to angle theta, up to precision 'prec', in polar coords?
  * Uses unit inner product with compact trig identities to reduce cost.  
  * Checks for azimuth alignment, which must be less than or equal to total precision, 
  * in order to prevent uncontrolled region growth near the vertical.
  * Note that lines which do have gradients along this direction may fragment, but
  * adjacent edges will grow as expected, mitigating overall fragmentation (unconfirmed). 
  */
int isaligned3(double grads_az,double grads_el,double theta_az,double theta_el,double prec)
{
  if( prec < 0.0 ) error("isaligned3: 'prec' must be positive.");
  if( grads_az == NOTDEF ) return FALSE;
  if( grads_el == NOTDEF ) return FALSE;  

  double sGel,sTel,cGel,cTel;
  sincos(grads_el,&sGel,&cGel);
  sincos(theta_el,&sTel,&cTel);
  
  /*Inner product, measuring parallel alignment*/
  double diff = sGel*sTel*cos(grads_az-theta_az) + cGel*cTel;
  //return fabs(diff)>=cos(prec);
  return diff>=cos(prec);
}

/*----------------------------------------------------------------------------*/
/** Is point (x,y,z) aligned to angle theta, up to precision 'prec', in polar coords?
  * If so, return sign of inner product of alignment, for parallel/anti-parallel analysis
  */
int isaligned3_sign(double grads_az,double grads_el,double theta_az,double theta_el,double prec)
{
  if( prec < 0.0 ) error("isaligned3: 'prec' must be positive.");
  if( grads_az == NOTDEF ) return FALSE;
  if( grads_el == NOTDEF ) return FALSE;  

  double sGel,sTel,cGel,cTel;
  sincos(grads_el,&sGel,&cGel);
  sincos(theta_el,&sTel,&cTel);
  
  /*Inner product, measuring parallel alignment*/
  double diff = sGel*sTel*cos(grads_az-theta_az) + cGel*cTel;
  return diff>=0;
}


/*----------------------------------------------------------------------------*/
/** Is point (x,y,z) aligned to angle theta, up to precision 'prec', in polar coords?
  * Uses unit inner product with compact trig identities to reduce cost.  
  * Checks for azimuth alignment, which must be less than or equal to total precision, 
  * in order to prevent uncontrolled region growth near the vertical.
  * Note that lines which do have gradients along this direction may fragment, but
  * adjacent edges will grow as expected, mitigating overall fragmentation (unconfirmed). 
  */
int isaligned3ORTH(double grads_az,double grads_el,double theta_az,double theta_el,double prec)
{
  if( prec < 0.0 ) error("isaligned3: 'prec' must be positive.");
  if( grads_az == NOTDEF ) return FALSE;
  if( grads_el == NOTDEF ) return FALSE;  
  
  double sGel,sTel,cGel,cTel;
  sincos(grads_el,&sGel,&cGel);
  sincos(theta_el,&sTel,&cTel);
  
  /*Inner product, measuring parallel alignment*/
  double diff = sGel*sTel*cos(grads_az-theta_az) + cGel*cTel;
  // azimithal alignment unrequired for othogonality measure
  return (1.-fabs(diff))>=cos(prec);
}


/*----------------------------------------------------------------------------*/
/** Switch orientation to antiparallel orientation 
  */
void align3(double * az,double *  el)

{
  (*az) -= M_PI;
  (*el) -= M_PI/2.;
  if ((*az) <= -M_PI) (*az) += M_PI*2.;
  if ((*el) < 0) (*el) += M_PI;
}
