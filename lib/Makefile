# MAKEFILE for the PyCIS-LSD analysis in C
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
# define flags (note CFLAGS -MMD and -I$(INCLUDE))
# some systems may need to specify 'gcc -std=c99'
CC = gcc 
INCLUDES = ../gsl/include  
LIBRARIES = ../gsl/lib 
CFLAGS = -MMD -fPIC -Wall -Wextra -O3 -g -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-unused-function -Wno-sign-compare -fopenmp -I$(INCLUDES)
LDFLAGS = -shared -fopenmp -lgsl -lgomp -lgslcblas -fopenmp
RM = rm -f
TARGET_LIB = libpycis_all.so

# list source files and generate objects and dependencies (via MMD)
SRCS = $(wildcard *.c) 
OBJS = $(SRCS:.c=.o)
DEPS = $(OBJS:.o=.d)

# define target shared library object 
# "make" is a default alias to "make all"
.PHONY: all
	all: %{TARGET_LIB}

# build target @ by applying ldflags with MMD to all objects
$(TARGET_LIB): $(OBJS)
	$(CC) -o $@ $^ ${LDFLAGS} -L$(LIBRARIES)

-include $(DEPS)

# clean up the temporary rule files and individual objects
# call this process by "make clean"
.PHONY: clean 
clean: 
	${RM} ${TARGET_LIB} ${OBJS}  ${DEPS}
