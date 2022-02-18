/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/markov.c: Determine first-order markov kernel for angular tolerance options

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
/*---------------------------------- Import ---------------------------------*/
/*----------------------------------------------------------------------------*/


#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include <tgmath.h>
#include <limits.h>
#include <float.h>
#include<string.h>
#include <time.h>
#include <gsl/gsl_sf_trig.h>
#include <sys/mman.h>

#include "markov.h"
#include "constants.h"
#include "misc.h"
#include "tuples.h"
#include "gradient.h"
#include "gaussians.h"

/*----------------------------------------------------------------------------*/
/*--------------------- Markov Transition Probabilities ----------------------*/
/*----------------------------------------------------------------------------*/

/* ----------------------------------------------------------------------------
  Compute markov transition matrix
  using "long double" data to mitigate round-off during evaluation 
*/
void NFA_matrix(double *output,double p0,double p11,double p01,int N)
{
  long double p1=  (long double) 1.0-p0;
  long double p10= (long double) 1.0-p11;
  long double p00= (long double) 1.0-p01;
    
  long double *plk0;
  plk0=(long double *) malloc(N*N*sizeof(long double));
  long double *plk1;
  plk1=(long double *) malloc(N*N*sizeof(long double));
  int i=0; int j=0;
  for(i=0;i<N;i++)
    {
      for(j=0;j<N;j++)
      {
        plk0[i*N+j]=(long double) 0;
        plk1[i*N+j]=(long double) 0;
        output[i*N+j]=(long double) 0;
      }
    }

  for(i=0;i<N;i++)
    for (j=0;j<N;j++)
      {
        if(i==0)
        {
          plk0[0+j]=(long double) 1;
          plk1[0+j]=(long double) 1;
        }
        else if(i==1)
        {
          plk0[i*N+j]=p01;
          plk1[i*N+j]=p11;
        }
        else
        {
          plk0[i*N+j]=(long double) 0;
          plk1[i*N+j]=(long double) 0;
        }
    }
  for(i=1;i<j;i++)
  {
    for(j=2;j<N;j++)
    {
      plk0[i*N+j]=plk0[i*N+j-1]*p00+plk1[(i-1)*N+j-1]*p01;
      plk1[i*N+j]=plk0[i*N+j-1]*p10+plk1[(i-1)*N+j-1]*p11;
    }
  }
  for(i=1;i<j;i++)
  {
    for(j=3;j<N;j++)
    {
      output[i*N+j]=(double) -log10l(plk0[i*N+j-1]*p0+plk1[(i-1)*N+j-1]*p1);
    } 
  }
  //free temp structures
  free((void *) plk0);
  free((void *) plk1);
}

/*----------------------------------------------------------------------------*/
/** Given a conditioning image, estimate the Markov transition probabilities  for 
 *  Gradient-by-Ratio computation, P(1|1) and P(1|0).  
 */
void make_markov( double * img, int X, int Y,
                           double ang_th, int n_bins,
                           double * inputv,double inputv_size)
{

  //fprintf(stdout,"MakeMarkov started\n");
  //fflush(stdout);

  //Instantiate variables per LSDSAR
  image_double image;
  image_double angles,modgrad;
  struct coorlist * list_p;
  void * mem_p;
  unsigned int xsize,ysize;
  double prec,p;
  /* check parameters */
  if( img == NULL || X <= 0 || Y <= 0 ) error("invalid image input.");
  if( ang_th <= 0.0 || ang_th >= 180.0 ) 
    error("'ang_th' value must be in the range (0,180).");
  if( n_bins <= 0 ) error("'n_bins' value must be positive.");
  /* angle tolerance */
  prec = M_PI * ang_th / 180.0;
  p = ang_th / 180.0;
  double beta;
  beta=inputv[0];

  /* load and scale image (if necessary) and compute angle at each pixel */
  image = new_image_double_ptr( (unsigned int) X, (unsigned int) Y, img );
  double scale,sigma_scale;
  scale=(double)inputv[5];
  sigma_scale=(double)inputv[6];
  if (scale != 1.0) 
  {
    image_double scaled_image;
    scaled_image = gaussian_sampler(image,scale,sigma_scale);
    angles = ll_angle( scaled_image, &list_p, &mem_p, &modgrad,
		       (unsigned int) n_bins,beta);
    free_image_double(scaled_image);
  }
  else angles = ll_angle( image, &list_p, &mem_p, &modgrad,
		       (unsigned int) n_bins,beta);
  xsize = angles->xsize;
  ysize = angles->ysize;
  
 //INDIVIDUAL LINE FREQUENCIES 
 int ang_mult, xx, yy;
 double hp_x0,hp_10,hp_x1,hp_11;
 double vp_x0,vp_10,vp_x1,vp_11;
 double hp_1,vp_1,hp_0,vp_0;
 int x_t, x_tminus;

 /* Calculate horizontal and vertical alignment freqencies*/
 // iterate precision options (prec, prec/2, prec/4)
 for(ang_mult=10;ang_mult<=14;ang_mult+=2)
 {
   hp_1=0; hp_0=0; vp_1=0;  vp_0=0; 
   // iterate vertical lines
   for(xx=0;xx<xsize;xx+=1)
   {
     vp_x0 = 0; vp_10 = 0; vp_x1 = 0; vp_11 = 0;
     for(yy=0;yy<ysize-1;yy+=1)
     {
        x_t        = isaligned(xx,  yy+1,  angles,0,prec);
        x_tminus   = isaligned(xx,  yy,    angles,0,prec); 
        if ( x_tminus==0 ) {++vp_x0; if ( x_t==1 ) ++vp_10;}
        else {++vp_x1; if ( x_t==1 ) ++vp_11;}
      }
      if( vp_x1>0) vp_1+=vp_11/vp_x1;//divide-by-zero risk 
      vp_0+=vp_10/vp_x0;
    }

   // iterate horizontal lines 
   for(yy=0;yy<ysize;yy+=1)
   {
     hp_x0 = 0; hp_10 = 0; hp_x1 = 0; hp_11 = 0;
     for(xx=0;xx<xsize-1;xx+=1)
     {
        x_t      = isaligned(xx+1,  yy,  angles,M_PI/2.,   prec);
        x_tminus = isaligned(xx,    yy,  angles,M_PI/2.,   prec);
        if ( x_tminus==0 ){++hp_x0; if ( x_t==1 ) ++hp_10;}
        else {++hp_x1; if ( x_t==1 ) ++hp_11;}              
      }
      if( hp_x1>0) hp_1+=hp_11/hp_x1;//divide-by-zero risk
      hp_0+=hp_10/hp_x0;
    }

     //Catch extrema cases 
     inputv[ang_mult]   = (hp_1 + vp_1)/(xsize+ysize);
     if(inputv[ang_mult]<=0) inputv[ang_mult]=0.0001;
     if(inputv[ang_mult]>=1) inputv[ang_mult]=0.9999;
     inputv[ang_mult+1] = (hp_0 + vp_0)/(xsize+ysize);
     if(inputv[ang_mult+1]<=0) inputv[ang_mult+1]=0.0001;
     if(inputv[ang_mult+1]>=1) inputv[ang_mult+1]=0.9999;
      
     // reduce tolerance for next loop 
     prec/=2;
     
 }
  
  free( (void *) image );   /* only the double_image structure should be freed,
                               the data pointer was provided to this functions
                              and should not be destroyed.                 */
  free_image_double(angles);
}


/*----------------------------------------------------------------------------*/
/** Given a conditioning image, estimate the Markov transition probabilities  for 
 *  Gradient-by-Ratio computation, P(1|1) and P(1|0).  
 */
void make_markovORTH( double * img, int X, int Y,
                           double ang_th, int n_bins,
                           double * inputv,double inputv_size)
{

  //fprintf(stdout,"MakeMarkov started\n");
  //fflush(stdout);

  //Instantiate variables per LSDSAR
  image_double image;
  image_double angles,modgrad;
  struct coorlist * list_p;
  void * mem_p;
  unsigned int xsize,ysize;
  double prec,p;
  /* check parameters */
  if( img == NULL || X <= 0 || Y <= 0 ) error("invalid image input.");
  if( ang_th <= 0.0 || ang_th >= 180.0 ) 
    error("'ang_th' value must be in the range (0,180).");
  if( n_bins <= 0 ) error("'n_bins' value must be positive.");
  /* angle tolerance */
  prec = M_PI * ang_th / 180.0;
  p = ang_th / 180.0;
  double beta;
  beta=inputv[0];

  /* load and scale image (if necessary) and compute angle at each pixel */
  image = new_image_double_ptr( (unsigned int) X, (unsigned int) Y, img );
  double scale,sigma_scale;
  scale=(double)inputv[5];
  sigma_scale=(double)inputv[6];
  if (scale != 1.0) 
  {
    image_double scaled_image;
    scaled_image = gaussian_sampler(image,scale,sigma_scale);
    angles = ll_angle( scaled_image, &list_p, &mem_p, &modgrad,
		       (unsigned int) n_bins,beta);
    free_image_double(scaled_image);
  }
  else angles = ll_angle( image, &list_p, &mem_p, &modgrad,
		       (unsigned int) n_bins,beta);
  xsize = angles->xsize;
  ysize = angles->ysize;
  
 //INDIVIDUAL LINE FREQUENCIES 
 int ang_mult, xx, yy;
 double hp_x0,hp_10,hp_x1,hp_11;
 double vp_x0,vp_10,vp_x1,vp_11;
 double hp_1,vp_1,hp_0,vp_0;
 int x_t, x_tminus;

 /* Calculate horizontal and vertical alignment freqencies*/
 // iterate precision options (prec, prec/2, prec/4)
 for(ang_mult=10;ang_mult<=14;ang_mult+=2)
 {
   hp_1=0; hp_0=0; vp_1=0;  vp_0=0; 
   // iterate vertical lines
   for(xx=0;xx<xsize;xx+=1)
   {
     vp_x0 = 0; vp_10 = 0; vp_x1 = 0; vp_11 = 0;
     for(yy=0;yy<ysize-1;yy+=1)
     {
        x_t        = isalignedORTH(xx,  yy+1,  angles,0,prec);
        x_tminus   = isalignedORTH(xx,  yy,    angles,0,prec); 
        if ( x_tminus==0 ) {++vp_x0; if ( x_t==1 ) ++vp_10;}
        else {++vp_x1; if ( x_t==1 ) ++vp_11;}
      }
      if( vp_x1>0) vp_1+=vp_11/vp_x1;//divide-by-zero risk 
      vp_0+=vp_10/vp_x0;
    }

   // iterate horizontal lines 
   for(yy=0;yy<ysize;yy+=1)
   {
     hp_x0 = 0; hp_10 = 0; hp_x1 = 0; hp_11 = 0;
     for(xx=0;xx<xsize-1;xx+=1)
     {
        x_t      = isalignedORTH(xx+1,  yy,  angles,M_PI/2.,   prec);
        x_tminus = isalignedORTH(xx,    yy,  angles,M_PI/2.,   prec);
        if ( x_tminus==0 ){++hp_x0; if ( x_t==1 ) ++hp_10;}
        else {++hp_x1; if ( x_t==1 ) ++hp_11;}              
      }
      if( hp_x1>0) hp_1+=hp_11/hp_x1;//divide-by-zero risk
      hp_0+=hp_10/hp_x0;
    }

     //Catch extrema cases 
     inputv[ang_mult]   = (hp_1 + vp_1)/(xsize+ysize);
     if(inputv[ang_mult]<=0) inputv[ang_mult]=0.0001;
     if(inputv[ang_mult]>=1) inputv[ang_mult]=0.9999;
     inputv[ang_mult+1] = (hp_0 + vp_0)/(xsize+ysize);
     if(inputv[ang_mult+1]<=0) inputv[ang_mult+1]=0.0001;
     if(inputv[ang_mult+1]>=1) inputv[ang_mult+1]=0.9999;
      
     // reduce tolerance for next loop 
     prec/=2;
     
 }
  
  free( (void *) image );   /* only the double_image structure should be freed,
                               the data pointer was provided to this functions
                              and should not be destroyed.                 */
  free_image_double(angles);
}


/*----------------------------------------------------------------------------*/
/*------------------- 3D-space Markov Transition Matrices --------------------*/
/*----------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
/*Fast orthogonality checks (reduced trig functions) */
static int isaligned3_markovVORTH(double grads_az,double grads_el,double cprec)
{
  if( grads_az == NOTDEF || grads_el == NOTDEF ) return FALSE;
  double diff = fabs(cos(grads_az)*sin(grads_el));
  return (1.-diff) >= cprec;
}
static int isaligned3_markovHORTH(double grads_az,double grads_el,double cprec)
{
  if( grads_az == NOTDEF || grads_el == NOTDEF ) return FALSE;
  double diff = fabs(sin(grads_az)*sin(grads_el))    ;
  return (1.-diff) >= cprec;
}
static int isaligned3_markovDORTH(double grads_az,double grads_el,double cprec)
{
  if( grads_az == NOTDEF || grads_el == NOTDEF ) return FALSE;
  double diff =  fabs(cos(grads_el));
  return (1.-diff) >= cprec;
}

/*------------------------------------------------------------------------*/
/*Fast alignment checks (reduced trig functions) */
static int isaligned3_markovV(double grads_az,double grads_el,double cprec)
{
  if( grads_az == NOTDEF || grads_el == NOTDEF ) return FALSE;
  double diff = (cos(grads_az)*sin(grads_el));
  //double diff = (cos(grads_az)*sin(grads_el));
  return (diff) >= cprec;
}
static int isaligned3_markovH(double grads_az,double grads_el,double cprec)
{
  if( grads_az == NOTDEF || grads_el == NOTDEF ) return FALSE;
  double diff = (sin(grads_az)*sin(grads_el))    ;
  //double diff = (sin(grads_az)*sin(grads_el))    ;
  return (diff) >= cprec;
}
static int isaligned3_markovD(double grads_az,double grads_el,double cprec)
{
  if( grads_az == NOTDEF || grads_el == NOTDEF ) return FALSE;
  double diff =  (cos(grads_el));
  //double diff =  (cos(grads_el));
  return (diff) >= cprec;
}

/*------------------------------------------------------------------------*/
/** Given a conditioning image, estimate the Markov transition probabilities  for 
 *  Gradient-by-Ratio computation, P(1|1) and P(1|0).  
 *
 *  MPI parallelization is used for acceleration 
 *
 */
void make_markov3( double * img, int X, int Y, int Z,
                          double ang_th, int n_bins, double * inputv,double inputv_size, 
                          int orth)
{

  //fprintf(stdout,"MakeMarkov3 started\n");
  //fflush(stdout);

  //Instantiate variables per LSDSAR
  image3_double image;
  image3_double modgrad;
  grads angles; 
  struct coorlist3 * list_p;
  void * mem_p;
  unsigned int xsize,ysize,zsize;
  double prec,p;
  /* check parameters */
  if( img == NULL || X <= 0 || Y <= 0 || Z<= 0 ) error("invalid image input.");
  if( ang_th <= 0.0 || ang_th >= 180.0 ) 
    error("'ang_th' value must be in the range (0,180).");
  if( n_bins <= 0 ) error("'n_bins' value must be positive.");
  if ((orth!=0) && (orth!=1)) error("make_markov3: post value not 0/1");
  
  /* angle tolerance */
  prec = M_PI * ang_th / 180.0;
  p = ang_th / 180.0;
  double beta;
  beta=inputv[0];
  if(prec<0.0) error("MakeMarkov3: 'prec' must be positive");

  /* load and scale image (if necessary) and compute angle at each pixel */
  image = new_image3_double_ptr((unsigned int)X, (unsigned int)Y, (unsigned int)Z, img);  
  double start,end;
  double scale,sigma_scale;
  scale=(double)inputv[5];
  sigma_scale=(double)inputv[6];
  
  if (scale != 1.0) 
  {
    image3_double scaled_image;
    scaled_image = gaussian3_sampler(image,scale,sigma_scale);
    angles = ll_angle3( scaled_image, &list_p, &mem_p, &modgrad,
		       (unsigned int) n_bins,beta);
    free_image3_double(scaled_image);
  }
  else angles = ll_angle3( image, &list_p, &mem_p, &modgrad,
		       (unsigned int) n_bins,beta);

  xsize = angles->az->xsize;
  ysize = angles->az->ysize;
  zsize = angles->az->zsize;

  //INDIVIDUAL LINE FREQUENCIES 
  //note that 'xx' is an index, while 'x' is the binary 0/1 alignment value for the kernel
  int ang_mult, xx, yy, zz;
  double hp_x0,hp_10,hp_x1,hp_11;
  double vp_x0,vp_10,vp_x1,vp_11;
  double dp_x0,dp_10,dp_x1,dp_11;
  double hp_1,vp_1,dp_1,hp_0,vp_0,dp_0;
  int x_t, x_tminus;

  angles3 angv,angh,angd;
  /* dl={cos(az)sin(el),sin(az)sin(el),cos(el)}
  * hence el=pi/2 for the horizontal case (sin=1,cos=0)
  * and   el=0    for the vertical   case (sin=0,cos=1)*/
  angv = new_angles3(0.,M_PI/2.);
  angh = new_angles3(M_PI/2.,M_PI/2.);//.);
  angd = new_angles3(0.,0.);
  /*temp data access copy*/
  double * azimg;
  azimg = (double *) malloc(xsize*ysize*zsize*sizeof(double));
  double * elimg; 
  elimg = (double *) malloc(xsize*ysize*zsize*sizeof(double));
  int i;
  for(i=0;i<(xsize*ysize*zsize);i++)
  {
    azimg[i]= angles->az->data[i];
    elimg[i]= angles->el->data[i];
  }
  double grads_az,grads_el;
  double cprec,cel,sel;
 
  /*Estimate markov kernel using OMP acceleration */
  int ic=0;
  if (inputv_size>16) 
  {
    ic=6;
    //printf("averaging 2D kernels...");fflush(stdout);
    double * imgsub = (double *) malloc(xsize*ysize*sizeof(double)); 
    double temp1,temp2,temp3,temp4,temp5,temp6;
    temp1=0; temp2=0; temp3=0; temp4=0; temp5=0; temp6=0;
    int i,j,k;
    for(k=0; k<zsize; k++)
    {
      //load frame into passable memory 
      for(i=0;i<xsize;i++)
        for(j=0;j<ysize;j++)
          imgsub[i+j*xsize] = img[k+zsize*(i+j*xsize)];
      if (orth==0)
        make_markov(imgsub, X, Y, ang_th, n_bins,inputv,inputv_size);
      else
        make_markovORTH(imgsub, X, Y, ang_th, n_bins,inputv,inputv_size);
      temp1+=(double)inputv[10];
      temp2+=(double)inputv[11];
      temp3+=(double)inputv[12];
      temp4+=(double)inputv[13];
      temp5+=(double)inputv[14];
      temp6+=(double)inputv[15];
    }
    free(imgsub);
    temp1/=zsize; temp2/=zsize; temp3/=zsize;
    temp4/=zsize; temp5/=zsize; temp6/=zsize;
    inputv[10]=temp1; inputv[11]=temp2; inputv[12]=temp3;
    inputv[13]=temp4; inputv[14]=temp5; inputv[15]=temp6;
  }
    
  for(ang_mult=10+ic;ang_mult<=14+ic;ang_mult+=2)
  {
    //Clear counting data  
    hp_1=0.; hp_0=0.; vp_1=0.;  vp_0=0.; dp_1=0.; dp_0=0.;
    cprec=cos(prec);
    // iterate vertical lines
    start=omp_get_wtime();
    #pragma omp parallel default(none) shared(xsize,ysize,zsize,azimg,elimg,prec,cprec,vp_0,vp_1,orth) private(xx,yy,zz,grads_az,grads_el,x_t,x_tminus,vp_x1,vp_x0,vp_11,vp_10)   
    {
    #pragma omp for reduction(+:vp_0) reduction(+:vp_1)
    for(xx=0;xx<xsize;xx+=1)
    {
      vp_x0 = 0; vp_10 = 0; vp_x1 = 0; vp_11 = 0;
      for(yy=0;yy<ysize-1;yy+=1)
      for(zz=0;zz<zsize-1;zz+=1)
      {      
        //next point in iteration should be along inner loop (z)    
        grads_az = azimg[ (zz+1) + zsize*((xx+0) + (yy+0) * xsize) ]; 
        grads_el = elimg[ (zz+1) + zsize*((xx+0) + (yy+0) * xsize) ]; 
        if(orth==0) x_t = isaligned3_markovV(grads_az,grads_el,cprec);
        else        x_t = isaligned3_markovVORTH(grads_az,grads_el,cprec);
        grads_az = azimg[ (zz) + zsize*((xx) + (yy) * xsize) ]; 
        grads_el = elimg[ (zz) + zsize*((xx) + (yy) * xsize) ]; 
        if(orth==0) x_tminus=isaligned3_markovV(grads_az,grads_el,cprec);                
        else        x_tminus=isaligned3_markovVORTH(grads_az,grads_el,cprec);
        if ( x_tminus==0 ) {++vp_x0; if ( x_t==1 ) ++vp_10;}
        else {++vp_x1; if ( x_t==1 ) ++vp_11;}
      } 
      if( vp_x1>0) vp_1+=vp_11/vp_x1;//divide-by-zero risk
      vp_0+=vp_10/vp_x0;
    }
    }
    end=omp_get_wtime();

    #pragma omp parallel default(none) shared(xsize,ysize,zsize,azimg,elimg,prec,cprec,hp_0,hp_1,orth) private(xx,yy,zz,grads_az,grads_el,x_t,x_tminus,hp_x1,hp_x0,hp_11,hp_10)   
    {
    #pragma omp for reduction(+:hp_0) reduction(+:hp_1)
    for(yy=0;yy<ysize;yy+=1)
    {
      hp_x0 = 0; hp_10 = 0; hp_x1 = 0; hp_11 = 0;
      for(xx=0;xx<xsize-1;xx+=1)
        for(zz=0;zz<zsize-1;zz+=1)
        {
          //next increment along z (inner loop)      
          grads_az = azimg[ (zz+1) + zsize*((xx+0) + (yy+0) * xsize) ]; 
          grads_el = elimg[ (zz+1) + zsize*((xx+0) + (yy+0) * xsize) ]; 
          if(orth==0) x_t = isaligned3_markovH(grads_az,grads_el,cprec);
          else        x_t = isaligned3_markovHORTH(grads_az,grads_el,cprec);
          grads_az = azimg[ (zz) + zsize*((xx) + (yy) * xsize) ]; 
          grads_el = elimg[ (zz) + zsize*((xx) + (yy) * xsize) ]; 
          if(orth==0) x_tminus = isaligned3_markovH(grads_az,grads_el,cprec);
          else        x_tminus = isaligned3_markovHORTH(grads_az,grads_el,cprec);
          if ( x_tminus==0 ){++hp_x0; if ( x_t==1 ) ++hp_10;}
          else {++hp_x1; if ( x_t==1 ) ++hp_11;}              
        }
      if( hp_x1>0) hp_1+=hp_11/hp_x1;
      hp_0+=hp_10/hp_x0;
     }
     }
     end=omp_get_wtime();
    
    #pragma omp parallel default(none) shared(xsize,ysize,zsize,azimg,elimg,cprec,dp_0,dp_1,orth) private(xx,yy,zz,grads_az,grads_el,x_t,x_tminus,dp_x1,dp_x0,dp_11,dp_10)   
    {
    #pragma omp for reduction(+:dp_0) reduction(+:dp_1)
    for(zz=0;zz<zsize;zz+=1)
    {
      dp_x0 = 0; dp_10 = 0; dp_x1 = 0; dp_11 = 0;
      for(xx=0;xx<xsize-1;xx+=1)
        for(yy=0;yy<ysize-1;yy+=1)
        {
          //NEXT increment along Y, (inner loop)
          grads_az=azimg[  (zz+0) + zsize*((xx+0) + (yy+1) * xsize)   ];
          grads_el=elimg[  (zz+0) + zsize*((xx+0) + (yy+1) * xsize)  ];
          if(orth==0) x_t = isaligned3_markovD(grads_az,grads_el,cprec);         
          else        x_t = isaligned3_markovDORTH(grads_az,grads_el,cprec);
          grads_az=azimg[  (zz) + zsize*((xx) + (yy) * xsize)   ];
          grads_el=elimg[  (zz) + zsize*((xx) + (yy) * xsize)  ];
          if(orth==0) x_tminus = isaligned3_markovD(grads_az,grads_el,cprec);
          else        x_tminus = isaligned3_markovDORTH(grads_az,grads_el,cprec);
          if ( x_tminus==0 ){++dp_x0; if ( x_t==1 ) ++dp_10;}
          else {++dp_x1; if ( x_t==1 ) ++dp_11;}              
        }
      if( dp_x1>0) dp_1+=dp_11/dp_x1;
      dp_0+=dp_10/dp_x0;
    }
    }
    end=omp_get_wtime();

    //Update inputv probability data
    int useXY=1; //use only xy or full xyz data for evaluation 
    if (inputv_size>16) useXY=0;
    //printf("vp (%.2f, %.2f) hp (%.2f, %.2f) dp (%.2f,%.2f)",vp_1,vp_0,hp_1,hp_0,dp_1,dp_0);fflush(stdout);
    if (useXY==1){
      hp_1 = (vp_1+hp_1)/(xsize+ysize); 
      hp_0 = (vp_0+hp_0)/(xsize+ysize);
    }
    else{
      if (useXY==-1){
        hp_1 = (dp_1)/(zsize); 
        hp_0 = (dp_0)/(zsize);
      }
      else {
        hp_1 = (vp_1+hp_1+dp_1)/(xsize+ysize+zsize); 
        hp_0 = (vp_0+hp_0+dp_0)/(xsize+ysize+zsize);
      }
    }
    inputv[ang_mult]   = hp_1;
    if(inputv[ang_mult]<=0) inputv[ang_mult]=0.0001;
    if(inputv[ang_mult]>=1) inputv[ang_mult]=0.9999;
    inputv[ang_mult+1] = hp_0;
    if(inputv[ang_mult+1]<=0) inputv[ang_mult+1]=0.0001;
    if(inputv[ang_mult+1]>=1) inputv[ang_mult+1]=0.9999;
    
    // reduce tolerance for next loop 
    prec/=2;   
  }
  free( (void *) image );   /* only the double_image structure should be freed,
                               the data pointer was provided to this functions
                              and should not be destroyed.                 */
  free_angles3(angh);
  free_angles3(angv);
  free_angles3(angd);
  free_grads(angles);
  //MAY REMOVE 
  free((void *)azimg);
  free((void *)elimg);
}
