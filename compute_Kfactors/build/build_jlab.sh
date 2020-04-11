#!/bin/bash

../configure CC="gcc -v" CXX="g++ -v" \
                CFLAGS="-O3  -fopenmp" CXXFLAGS="-std=c++0x -O3  -fopenmp" \
                --prefix="${HOME}/git/kinematic_factors_install" \
                --with-adat="${HOME}/git/adat_install" \
                --with-semble="${HOME}/git/semble_install" \
                #--with-gcc-8="/usr/local/Cellar/gcc/8.2.0/"
