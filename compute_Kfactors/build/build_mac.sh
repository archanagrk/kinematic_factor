#!/bin/bash

../configure --prefix="${HOME}/LQCDSoftware/Kinematic_factor/compute_Kfactors_install" \
             --with-adat="${HOME}/LQCDSoftware/Kinematic_factor/adat_devel_install" \
             --with-semble="${HOME}/LQCDSoftware/three_pt_analysis/semble_install" \
             CXXFLAGS="-O3 -mtune=native -fopenmp -fbounds-check" CFLAGS="-O3 -mtune=native -fopenmp -fbounds-check" \
             CXX=g++-9 CC=gcc-9
