#/bin/bash

cd build
rm -rf CMakeCache.txt  CMakeFiles  cmake_install.cmake  Makefile
cd ../bin 
rm -rf main
cd ..
cd build 
cmake ..
make
cd ..