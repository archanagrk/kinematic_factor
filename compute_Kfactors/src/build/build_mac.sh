#!/bin/bash

../configure --prefix="${HOME}/Kinematic_factor/compute_Kfactors_install" \
             --with-adat="${HOME}/Kinematic_factor/adat_devel_install" \
             CXXFLAGS="-O3 -mtune=native -fopenmp -fbounds-check" CFLAGS="-O3 -mtune=native -fopenmp -fbounds-check" \
             CXX=g++-8 CC=gcc-8
