#!/bin/bash
icc -c -O2 random.c
icc -c -O2 -std=c++11 MemTimeTester.cpp
icc -c -O2 -std=c++11 Matrix.cpp
icc -c -O2 -std=c++11 class_mc_io.cpp
icc -c -O2 -std=c++11 IsingLattice2D.cpp
icc -c -O2 -std=c++11 LongRangeWolff2D.cpp
nvcc -c -O2 -std=c++11 obs_calc_fast.cu
mpic++ -c -O2 -std=c++11 spin_boson_debug.cpp
mpic++ -std=c++11 -lcufft -lgsl -lgslcblas -I$GSL_INC -o spin_boson_debug.exe spin_boson_debug.o random.o class_mc_io.o LongRangeWolff2D.o IsingLattice2D.o Matrix.o MemTimeTester.o obs_calc_fast.o
