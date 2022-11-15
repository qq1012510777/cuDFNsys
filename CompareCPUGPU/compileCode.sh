#/bin/bash

if [[ -z "$1" ]]
then
    echo 'Please specify which cu file you want to compile (you can input 1, 2, 3, 4, 5 or 6)'
    exit 0
fi

cd build
# rm -rf CMakeCache.txt  CMakeFiles  cmake_install.cmake  Makefile
cd ../bin 
rm -rf main
cd ..
cd build 
cmake .. -DCUFILE=$1
make
cd ..