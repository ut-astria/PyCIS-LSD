# PyCIS-LSD : Python Computational Inference from Structure, Line Segment Detector


An a-contrario detection sub-algorithm for extracting narrow lines within dense optical data cubes.


Benjamin Feuge-Miller: benjamin.g.miller@utexas.edu
The University of Texas at Austin, 
Oden Institute Computational Astronautical Sciences and Technologies (CAST) group
*Date of Modification: December 23, 2021*

**NOTICE**: For copyright and licensing, see 'notices' at bottom of readme

------------------------------------------------------------------

## OVERVIEW:

This algorithm takes in a dense data cube of optical time-series data and looks for narrow lines corresponding to streaks or actively tracked trajectories.  
For a more detailed discussion of *a-contrario* analysis see the main PyCIS pipeline which calls this algorithm.


------------------------------------------------------------------

## OPERATION:

This software is designed to be used as a subsegment of the main PyCIS processing pipeline. 
The PyCIS main pipeline will control setup, if used.  Otherwise, see the main pipleine for instruction on operating PyCIS-LSD

**Setup**
   If scripts cannot be read, run:
   `<dos2unix *.sh>`
   Finally, install prerequisites and build the software by running:
   `<. setup.sh>`.



**Demo**:
  See the main PyCIS pipeline for demo scripts and data download instructions.  


**Output Files**:
Input and output locations and names are specified within the main file, see `<demo.py>` for example.

* Input: data/...
    * yyymmdd_norad_satname/*.fit - a folder with raw fits frames [see citations in DATASET section below] 
* Output: 
    * results_work/... 
        * data1_x_y_name.npy - edge line detections, x and y labels for parallelization are absent when merged.
        * data2_x_y_name.npy - center line detections, x and y labels for parallelization are absent when merged.

Runtime estimates for the two-dataset demo are around 30 minutes on one TACC none with full memory, 48-thread opm capability.

------------------------------------------------------------------

## DEMO VISUAL OBSERVATION

See the main PyCIS pipeline for details.

------------------------------------------------------------------

## TOC:

* pycis_lsd_wrapper.py - interface for launching pycis.c through python
* pycis.c -              main python extension module 
* setup.py -             link python-c extension module
* setup.sh -             install gsl/pyenv/Makefile, link


* LIB:
    * constants -     set up some constants
    * gaussians -     gaussian down-sampling for antialiasing (unused)
    * gradient -      compute gradient magnitude/orientation, and alignment checks
    * interface -     pipeline helpers for pycis.c
    * lsd(2/3/3b) -   main lsd pipelines (spatial, temporal, and unified)
    * markov -        estimate markov kernels and build transition matrices
    * misc -          lists, angle functions
    * NFA -           estimate markov or negative binomial approximation tail probability
    * rectanges(2/3)- build rectangle objects and iterator to call nfa
    * regions(2/3) -  build and improve pixel regions for estimating rectangles
    * tuples -        construction of line tuples and image structures 

------------------------------------------------------------------

## NOTICES:


**NOTICE**: 
PyCIS-LSD: An a-contrario detection sub-algorithm for extracting narrow lines within dense optical data cubes.
Copyright (C) 2022, Benjamin G. Feuge-Miller, <benjamin.g.miller@utexas.edu>

PyCIS-LSD is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published 
by the Free Software Foundation, either version 3 of the License, 
or (at your option) any later version.

PyCIS-LSD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.
 
You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.


**NOTICE**: 
PyCIS-LSD is modified from the source code of LSDSAR:
"LSDSAR, a Markovian a contrario framework for line segment detection in SAR images"
by Chenguang Liu, Rémy Abergel, Yann Gousseau and Florence Tupin. 
(Pattern Recognition, 2019).
https://doi.org/10.1016/j.patcog.2019.107034
*Date of Modification: April 30, 2021*


**NOTICE**: 
LSDSAR is modified from the source code of LSD:
"LSD: a Line Segment Detector" by Rafael Grompone von Gioi,
Jeremie Jakubowicz, Jean-Michel Morel, and Gregory Randall,
Image Processing On Line, 2012. DOI:10.5201/ipol.2012.gjmr-lsd
http://dx.doi.org/10.5201/ipol.2012.gjmr-lsd
*Date of Modification: 27/06/2018*





------------------------------------------------------------------
