/*----------------------------------------------------------------------------
PyCIS - Python Computational Inference from Structure

lib/nfa.c: compute the NFA equation 
  using precomputed markov table and binomial approximation 

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
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_qrng.h>
#include<gsl/gsl_sf_gamma.h>
#include<gsl/gsl_sf_trig.h>
#include<sys/mman.h>

#include "nfa.h"
#include "constants.h"
#include "misc.h"
#include "tuples.h"

/*----------------------------------------------------------------------------*/
/*----------------------------- NFA computation ------------------------------*/
/*----------------------------------------------------------------------------*/


/*NOTE: Here cut vonGioi's log functions, unused in LSDSAR in favor of Markov table*/

/*----------------------------------------------------------------------------*/
/** Computes -log10(NFA).

    NFA stands for Number of False Alarms:
    @f[
        \mathrm{NFA} = NT \cdot B(n,k,p)
    @f]

    - NT       - number of tests
    - B(n,k,p) - tail of binomial distribution with parameters n,k and p:
    @f[
        B(n,k,p) = \sum_{j=k}^n
                   \left(\begin{array}{c}n\\j\end{array}\right)
                   p^{j} (1-p)^{n-j}% alpha=2;
*/
double nfa(int n, int k, double p, double logNT,double *mnfa,int N)
{
   if(n>N||k>N)
        return 101;
  /* check parameters */
  if( n<0 || k<0 || k>n || p<=0.0 || p>=1.0 )
    error("nfa: wrong n, k or p values.");
  /* trivial cases */
  if( n<3 || k==0 ) return -logNT;
  return mnfa[k*N+n]-logNT;
}

/*----------------------------------------------------------------------------*/
/** Computes -log10(NFA).
 * 
 *  Use a negative binomial approximation as formulated by 
 *  Xia, A. & Zhang, M. (2010). "On approximation of Markov binomial distributions".  
 *    Bernoulli 15(4), 2009, 1335-1350: DOI: 10.3150/09-BEJ194
*/
double nfaORTH(int n, int k, double pp, double logNT, double *mnfa, int N, double * divout)
{ 

  //let *mfa store :
  //     mnfa = [p0,p11,p01] for p for mfa, mnfa_2, mnfa_4     
  int nscale = 1;
  //if((n>N)||(k>N))
  //     return 101;
  /* check parameters */
  if( n<0 || k<0 || k>n || pp<=0.0 || pp>=1.0 )
    error("nfa: wrong n, k or p values.");
  /* trivial cases */
  if( n<3 || k==0 ) return -logNT;

  double p0  = mnfa[0];
  double p11 = mnfa[1];
  double p01 = mnfa[2];
  double p1=1.0-p0;
  double p10=1.0-p11;
  double p00=1.0-p01;
  double *plk0;

  double alpha = p01; 
  double beta = p11;
  double p  = alpha/(1.-beta+alpha);
  double pi = (1.-beta)/(1.-beta+alpha);

  double A0 = (2.*alpha*(1.-beta)*(beta-alpha))/pow(1.-beta+alpha,3.);
  double tempp2=pow(p,2.);
  double A1 = (2.*alpha*(1.-beta)*(beta-alpha))/pow(1.-beta+alpha,4.);
  double nn=(double) n;
  double ES   = nn*p;
  double VarS = nn*p*(1.-p)+nn*A0-A1+A1*pow(beta-alpha,nn);
  
  double output, div;
  double r, q, mhat, m, theta;
  
  double mu1 = (1.-alpha)/alpha;
  double sig1= (1.-alpha)/pow(alpha,2.);
  double mu2 = beta/(1.-beta);
  double sig2= beta/pow(1.-beta,2.);
  double C0 = fabs(beta-alpha)*(5.+max2(43.*alpha,beta))/pow(1.-max2(beta,alpha),2.);
  double C1 = 10.*max2(beta,alpha)/(1.-max2(beta,alpha));
  double C2 = (1.-p)*(5.+max2(23.*alpha,beta))/pow(1.-max2(alpha,beta),2.);
  double K1 = sqrt(5.)*sqrt((mu1+mu2+2.)/min2(min2(1.-alpha,beta),0.5));
  double K2 = 90.*(sig1+sig2)/(mu1+mu2+2.);
  
  double tempoutput=0;
  double temppost=0;
  if (VarS>=ES)  
  {
    r=pow(ES,2.) / (VarS-ES) ; 
    q=ES/VarS;
    temppost=pow(beta,floor(nn/4.));
    div = C0*(2.*K1/sqrt(nn)+4.*K2/nn + temppost);//pow(beta,floor(nn/4.)));
    /*output = 1 - NB(r,q) //negative binomial 
      gsl: (ushort k, p, n )
      match: p^n (1-p)^k = q^r (1-q)^k ; so (k,p,n)=(k,q,r) */
    tempoutput = gsl_cdf_negative_binomial_Q((unsigned int)k,(double)q,(double)r);//ushort k | p, n
    tempoutput=-log10(tempoutput)-logNT;
    if(!isfinite(tempoutput))
    {
      double to1a,to1b,to1c,to1,to2,to3;
      to1a= gsl_sf_lngamma((double)r+(double)k);
      to1b= -gsl_sf_lngamma((double)k+1.);//log10l(tgammal(kl+1.));
      to1c= -gsl_sf_lngamma((double)r);//log10l(tgammal(nl));
      to1 = to1a+to1b+to1c;
      to2 = log(q)*r ;
      to3 = log(1.-q)*(double)k; 
      tempoutput=(to1+to2+to3);// @p
      double tempoutput2;
      for(int kk=k+1;kk<floor(r);kk++)
      {
        //printf("TEST k %d, kk %d, floorr %d\n",(int)k,(int)kk,(int)floor(r));fflush(stdout);
        to1a= gsl_sf_lngamma((double)r+(double)kk);
        to1b= -gsl_sf_lngamma((double)kk+1.);//log10l(tgammal(kl+1.));
        to1c= -gsl_sf_lngamma((double)r);//log10l(tgammal(nl));
        to1 = to1a+to1b+to1c;
        //CONST:to2=-log10(q)*r ;
        to3=log(1.-q)*(double)kk; 
        tempoutput2=(to1+to2+to3);// @p
        //both temps are in -lot10 currently, convert back to basen before operation 
        tempoutput+=log1p(exp(tempoutput2-tempoutput));
      }  
      tempoutput/=-log(10);
      tempoutput-=logNT;
    }
  }
  else
  {
    mhat = pow(ES,2.)/(ES-VarS);
    m = floor(mhat);
    theta = nn*p/m;
    div = (C1*fabs(p-theta)/(1.-theta) + C2*fabs(beta-alpha)/(1.-theta))*(  2.*K1/sqrt(nn)+4.*K2/nn + pow(max2(beta,alpha),floor(n/4.))   )+(pow(theta,2.)*(mhat-m)/(nn*p*(1.-theta)));
    tempoutput = gsl_cdf_binomial_Q((unsigned int)k,(double)theta,(unsigned int)m);//ushort k | p, n
    tempoutput=-log10(tempoutput)-logNT;
    if(!isfinite(tempoutput)){// tempoutput=10001;
      double to1a,to1b,to1c,to1,to2,to3;
      to1a= gsl_sf_lngamma((double)m+1.);
      to1b= -gsl_sf_lngamma((double)k+1.);//log10l(tgammal(kl+1.));
      to1c= -gsl_sf_lngamma((double)m-(double)k+1.);//log10l(tgammal(nl));
      to1 = to1a+to1b+to1c;
      to2 = log(theta)*(double)k ;
      to3 = log(1.-theta)*((double)m-(double)k); 
      tempoutput=(to1+to2+to3);// @p
      double tempoutput2;
      for(int kk=k+1;kk<floor(m);kk++)
      {
        to1a= gsl_sf_lngamma((double)m+1.);
        to1b= -gsl_sf_lngamma((double)kk+1.);//log10l(tgammal(kl+1.));
        to1c= -gsl_sf_lngamma((double)m-(double)kk+1.);//log10l(tgammal(nl));
        to1 = to1a+to1b+to1c;
        to2 = log(theta)*(double)kk ;
        to3 = log(1.-theta)*((double)m-(double)kk); 
        tempoutput2=(to1+to2+to3);// @p
        //both temps are in -lot10 currently, convert back to basen before operation 
        tempoutput+=log1p(exp(tempoutput2-tempoutput));
      }  
      tempoutput/=-log(10);
      tempoutput-=logNT;
    }
  }
  if(!isfinite(tempoutput)) tempoutput=10001.;
  output = (double) tempoutput;
  *divout = div;
  return output;
}

