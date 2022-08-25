#/bin/bash

cd build
# rm -rf CMakeCache.txt  CMakeFiles  cmake_install.cmake  Makefile
cd ../bin 
rm -rf Transform2DH5ParticleDataTo3D
cd ..
cd build 
cmake ..
make
cd ..