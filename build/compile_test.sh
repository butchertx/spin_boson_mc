#!/bin/bash
gcc -c -O2 $SBLIB/random.c
g++ -I$SBLIB -c -O2 -std=c++11 $SBLIB/MemTimeTester.cpp
g++ -I$SBLIB -c -O2 -std=c++11 $SBLIB/Matrix.cpp
g++ -I$SBLIB -c -O2 -std=c++11 $SBLIB/class_mc_io.cpp
g++ -I$SBLIB -c -O2 -std=c++11 $SBLIB/IsingLattice2D.cpp
g++ -I$SBLIB -c -O2 -std=c++11 $SBLIB/LongRangeWolff2D.cpp
nvcc -I$SBLIB -c -O3 -std=c++11 $SBLIB/obs_calc_fast.cu
mpic++ -I$SBLIB -I$CPATH -c -O2 -std=c++11 $SBLIB/spin_boson_debug.cpp
mpic++ -std=c++11 -L$LD_LIBRARY_PATH -lcudart -lcufft -lgsl -lgslcblas -o ~/spin_boson_mc/bin/spin_boson_debug.exe spin_boson_debug.o random.o class_mc_io.o LongRangeWolff2D.o IsingLattice2D.o Matrix.o MemTimeTester.o obs_calc_fast.o
