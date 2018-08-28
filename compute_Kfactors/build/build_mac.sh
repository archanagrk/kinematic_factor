#!/bin/bash

../configure --prefix="${HOME}/Subduction/subduction_install" \
             --with-adat="${HOME}/Subduction/adat_install" \
             CXXFLAGS="-O3 -mtune=native -fopenmp -fbounds-check" CFLAGS="-O3 -mtune=native -fopenmp -fbounds-check" \
             CXX=g++-8 CC=gcc-8
