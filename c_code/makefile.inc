filein  = -I/usr/include -I/public/home/liuchang/src/testpycbc/opt/lscsoft/non-lsc/include

LDLIBS  = -L/usr/local/lib/

GSLIB =  -I/usr/local/include/ -lgsl -lgslcblas -lm

CXXAPPFLAGS = -O2 -Wno-deprecated -Dfortran3 -Dnewc -std=gnu++0x
f90appflags = -O2 -fpp
f90 = gfortran
f77 = g77
CXX = g++
CC  = gcc
CLINKER = g++
