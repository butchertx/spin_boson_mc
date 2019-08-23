#!/bin/bash
icc -c -O3 $SBLIB/random.c
icc -c -O3 -std=c++11 $SBLIB/Matrix.cpp
icc -c -O3 -std=c++11 $SBLIB/class_mc_io.cpp
icc -c -O3 -std=c++11 $SBLIB/LongRangeWolff2D.cpp
nvcc -c -O3 -std=c++11 $SBLIB/obs_calc_fast.cu
mpic++ -c -O3 -std=c++11 $SBLIB/spin_boson_mpi.cpp
mpic++ -std=c++11 -L$LD_LIBRARY_PATH -lcudart -lcufft -lgsl -lgslcblas -I$CPATH spin_boson_mpi.o random.o class_mc_io.o LongRangeWolff2D.o Matrix.o obs_calc_fast.o -o ../bin/spin_boson_mpi.exe 
