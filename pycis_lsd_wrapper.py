'''
PyCIS - Python Computational Inference from Structure

pylib/pycis_lsd_wrapper.py: Wraps the LSD functionality of pycis.c and lib/ files, 
    and provided switches for directing 2D/3D and centerline/edgeline calls. 

Benjamin Feuge-Miller: benjamin.g.miller@utexas.edu
The University of Texas at Austin, 
Oden Institute Computational Astronautical Sciences and Technologies (CAST) group
*Date of Modification: December 27, 2021

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
'''

## IMPORT NECESSARY LIBRARIES
import os
import numpy as np
import pylsd.pycis as pycis
import time
 
def get_size(I):
    '''    Get matrix size variables with necesary 0-setting    '''
    M = len(I)
    if (I.ndim==1):
        M=len(I)
        N=0
        O=0
    elif (I.ndim==2) and (M>10):
        M,N=np.shape(I)
        O=0
    elif (I.ndim==3) and (M>10):
        M,N,O=np.shape(I)
    else:
        M,N,O=0,0,0
    return M,N,O

def get_sizes(I,I0):
    '''    Get matrix size variables with necessary 0-settings    '''
    M,N,O = get_size(I)
    #print('M %d, N %d, O %d'%(M,N,O))
    M0,N0,O0 = get_size(I0)
    #print('M0 %d, N0 %d, O0 %d'%(M0,N0,O0))
    return M,N,O,M0,N0,O0

def my_unflatten(I,M,N,O):
    '''    Format image data for passing to pycis    '''
    if (O==0):
        templist = np.zeros((M,N))
        for i in range(M):
            for j in range(N):
                templist[i,j]=I[(i+M*j)]
    if (O>0):
        templist = np.zeros((M,N,O))
        for i in range(M):
            for j in range(N):
                for k in range(O):
                    templist[i,j,k]=I[k+O*(i+M*j)]
        I = templist
    return I

def my_flatten(I,M,N,O):
    '''    Format image data for passing to pycis    '''
    if (O==0):
        I = I.T.flatten().tolist()
    if (O>0):
        templist = np.zeros((M*N*O,))
        for i in range(M):
            for j in range(N):
                for k in range(O):
                    templist[k+O*(i+M*j)] = I[i,j,k]
        I = templist.tolist()
    return I


def main(I,I0, folder='results',name='temp',
    a=4.,d=.4,t=1.,scale=0.8,sigma=0.6,e=0,getp=1,p=[0,0,0,0,0,0,0,0,0,0,0,0],p2=[0,0,0,0,0,0,0,0,0,0,0,0],idxs=[-1,-1],shape=[4096,4096,25]):
    '''
    Run LSD pipeline for a test and conditioning image pair .
    Reshape data and set input parameter vector.
    Decide between output functions.
    
    Input:  
        I :     test image
        I0:     noise model (or prior edge lines)
        folder: location for results 
        name:   name for saving .png and .npy files
        a:      gradient-by-ratio parameter.  Produces a (k*2+1)^d kernel for k=log(10)*a
                    eg a=1 produces a 7^d kernel, a=2 an 11^d kernel, a=3 a 15^d kernel, etc
        d:      density threshold for improving regions.  Higher values enforce stronger linearity constraints   
        t:      denominator factor for tolerance.  Initial tolerance is 22.5deg, empirically optimal for 2D
        p:      markov kernel for estimating parallel alignments 
                    (in edge case, build regions by parallel alignment and count alignments )
        p2:     markov kernel for estimating parallel alignments 
                    (in edge case, build regions by parallel alignment and count alignments )
        getp:   flag for pipline - control output as markov kernel or lines.  see pycis.c
        scale:  fraction of gaussian downsampled side lengths to input data (set to 1 for no downsampling)
        sigma:  std=sigma/scale for variance of gaussian downsampling via seperable 1D kernels
        e:      -log10(epsilon) NFA threshold.  Choose 0 for epsilon=1 (default theory).  
                    Robust to selection, but setting very large (e=6) can reduce spurious noise
                    and improve 2nd-order detection if there exists a statistically relevant number of detections
           
    Output: one of the following: 
        inputv/inputvorth:  updated settings vector with solved markov kernels
        data1_name.npy:     edge line detections
        data2_name.npy:     center line detections

    Notes: I,I0 cannot be empty vectors.  To set an 'empty image',
           use a small matrix, e.g. 2x2 identity, and the function
           will set the variables X,X0 appropriatly.
    '''
    if not os.path.exists(folder):
        os.makedirs(folder)
    '''
    #global I3
    #global I03
    x1=I[0]; x2=I[1]; y1=I[2]; y2=I[3]
    try:
        I = I3[x1:x2, y1:y2, :]
    except: 
        I = I3[x1:x2, :]
    x1=I[0]; x2=I[1]; y1=I[2]; y2=I[3]
    try:
        I0 = I03[x1:x2, y1:y2, :]
    except: 
        I0 = I03[x1:x2,:]
    '''
    ## SET INPUT DIRECTIONS FOR VERSIONING - LSD CANNOT ACCEPT EMPTY IMAGES, BUT REQUIRES M/M0=0
    M,N,O,M0,N0,O0 = get_sizes(I,I0)
  
    ## FLATTEN IMAGES FOR USE IN C, COPY FOR PLOTTING 
    I_full = np.copy(I)
    I0_full= np.copy(I0)
    I = my_flatten(I,M,N,O)
    I0= my_flatten(I0,M0,N0,O0)
    
    ## SET INPUT PARAMETER VECTOR 
    p11= p[0];p01= p[1];p11_2= p[2];p01_2= p[3];p11_4= p[4];p01_4= p[5]
    dp11=p2[0];dp01=p2[1];dp11_2=p2[2];dp01_2=p2[3];dp11_4=p2[4];dp01_4=p2[5]
    pb11= p[6];pb01= p[7];pb11_2= p[8];pb01_2= p[9];pb11_4= p[10];pb01_4= p[11]
    dpb11=p2[6];dpb01=p2[7];dpb11_2=p2[8];dpb01_2=p2[9];dpb11_4=p2[10];dpb01_4=p2[11]
    alpha=a#4.
    eps=10.**(-1.*e)#(1/1.)
    density=d #0.4
    angth=22.5/t
    sizenum=np.sqrt(M**2.+N**2.)*5.
    if sizenum>(10.**4):
        sizenum=10.**4
    if O>0:
        sizenum=min(10.**6.,np.sqrt(M**2.+N**2.+O**2.)*5.)
    inputv=[alpha,eps,density,sizenum,angth,
        scale,sigma,shape[0],shape[1],shape[2],
        p11, p01, p11_2, p01_2, p11_4, p01_4,
        pb11, pb01, pb11_2, pb01_2, pb11_4, pb01_4]
    inputvorth=[alpha,eps,density,sizenum,angth,
        scale,sigma,shape[0],shape[1],shape[2],
        dp11,dp01,dp11_2,dp01_2,dp11_4,dp01_4,
        dpb11,dpb01,dpb11_2,dpb01_2,dpb11_4,dpb01_4]
   
    ## RUN LSD
    #Get markov kernel if requested
    markov = getp
    if markov>1:
        if markov==2:
            #print('\n--------------- MARKOV-PARALLEL: %s ---------------\n'%name,flush=True)              
            print('MARKOV-PARALLEL: %s'%name,flush=True)  
        else:
            #print('\n--------------- MARKOV-ORTHOGONAL: %s ---------------\n'%name,flush=True)  
            print('MARKOV-ORTHOGONAL: %s'%name,flush=True)  

        time.sleep(1)
        lines = pycis.pycis(I,M,N,O,I0,M0,N0,O0,inputv,inputvorth,markov)
        del I, I0, alpha, density, angth, a,d,t,inputv,I_full,I0_full
        #Return inputv 
        if (idxs[0]==-1) and (idxs[1]==-1):
            return lines 
        else:
            return (idxs,lines)   
    else:
        #Find kernel and run estimation...
        savename='%s/data1_%s'%(folder,name)
        if markov==1: #markov = 1, run markov estimation plus edge detection 
            #print('\n--------------- LSD + Markov: %s ---------------\n'%name,flush=True)    
            print('LSD + Markov: %s'%name,flush=True)    
        elif markov==0 and M0==0: #markov=0, edges with prior markov
            #print('\n--------------- LSD (Edge): %s ---------------\n'%name,flush=True)      
            print('LSD (Edge): %s'%name,flush=True)      
        elif markov==0 and M0>0: #markov=0, centerlines with prior markov
            #print('\n--------------- LSD (Centers): %s ---------------\n'%name,flush=True)   
            print('LSD (Centers): %s'%name,flush=True)   
            savename = '%s/data2_%s'%(folder,name) 
        elif markov<0: #markov = 1, run markov estimation plus edge detection 
            pass
        else:
            print("ERROR: incompatible markov and X0 input .",flush=True)
            quit()

        time.sleep(1)
        lines = pycis.pycis(I,M,N,O,I0,M0,N0,O0,inputv,inputvorth,markov)
        if markov<0: #markov = 1, run markov estimation plus edge detection 
            return my_unflatten(np.asarray(lines),M,N,O)
        #lines = [];

        ## SAVE LINE RESULTS 
        np.save(savename,lines)
      
        del I, I0, alpha,density,angth, a,d,t,inputv,I_full,I0_full
        #Return nothing - data is saved to file
        return lines