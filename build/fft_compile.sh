#!/bin/bash

CPPFLAGS="-std=c++11 -DTEST=1"

icc $CPPFLAGS -o fft_test.x ../lib/fftw_helper.cpp -lfftw3 -lm
