"""
Setup function to link pycis.c, libpycis_all.so, and gsl 
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
"""

from distutils.core import setup, Extension 
import os
import sysconfig

d = os.getcwd()
exargs = sysconfig.get_config_var('CFLAGS').split()
exargsL = sysconfig.get_config_var('LDFLAGS').split()
#Some systems may need to add '-std=c99'
exargs += ['-fopenmp','-O3',
           "-Wno-unused-variable","-Wno-unused-but-set-variable",
           "-Wno-unused-function","-Wno-sign-compare"]
exargsL += ['-fopenmp']
#gsllib = os.environ.get('TACC_GSL_LIB') 
#gslinc = os.environ.get('TACC_GSL_INC')
gsllib ='%s/gsl/lib'%d 
gslinc ='%s/gsl/include'%d 
module = Extension('pycis',
                    include_dirs = [gslinc, 
                                    '%s/lib'%d,
                                    'gen','two','three'],
                    libraries    = ['gsl', 'gslcblas','gomp', 'pycis_all'],
                    library_dirs = [gsllib,
                                    '%s/lib'%d,],
                    runtime_library_dirs = [gsllib,
                                    '%s/lib'%d],
                    sources = ['pycis.c'],
                    extra_compile_args = exargs,
                    extra_link_args    = exargsL)

setup(name = "pycis", ext_modules = [module])
