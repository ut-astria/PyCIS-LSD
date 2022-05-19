#!/bin/bash
#
#This sofware should be located within a sub-folder of the main pycis directory 
#
#Download, install, clean up, and set paths for:
#   1) GSL 2.6
#   2) Python3 environment (append to pycis environment one level higher)
#Install astrometry index files for solving
#Then make the C library and launch python setup to link PyCIS
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


if type "sbatch" >& /dev/null; then
    module load gcc/9.1
    module load cmake/3.16.1 #10.2
    module load python3/3.8.2
fi

export CC=`which gcc`
export CFLAGS="-O3 -march=native"

CWDLSD=$(pwd)
#export EXTRA_CFLAGS=" -std=c99 "
# install gsl
if [ ! -d "./gsl" ]; then 
    # get gsl
    wget -c ftp://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz 
    tar -zxf gsl-2.6.tar.gz
    rm gsl-2.6.tar.gz
    # prepare paths
    mkdir gsl
    cd gsl-2.6
    # install gsl
    ./configure --prefix=${CWDLSD}/gsl
    if type "sbatch" >& /dev/null; then
        make
    else
        make -j 12
    fi
    make check
    make install
    # clean up
    cd ${CWDLSD}
    rm -rf gsl-2.6
    
fi
export PATH=$PATH:$CWDLSD/gsl/bin
export PATH=$PATH:$CWDLSD/gsl/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CWDLSD/gsl/lib
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:$CWDLSD/gsl/lib/pkgconfig/

# install venv using gsl, to install astropy
#some systems may need python3.6 specifically if 3.8+ is not available 
# USE ENV FROM ONE DIRECTORY UP
if [ ! -d "../env" ]; then
    python3 -m venv ../env
fi
if [ ! -d "../env/lib/python3.8/site-packages/numpy" ]; then
    source ../env/bin/activate
    python3 -m pip install --upgrade pip
    python3 -m pip install -r requirementslsd.txt
    deactivate
fi

#activate env
source ../env/bin/activate

#PyCIS: Install  Astrometry.net, and CFITSIO/WCLIB

#install library sources
cd lib
make
cd ${CWDLSD}

#python3 setup.py install --user
python3 setuplsd.py build_ext --inplace

deactivate

