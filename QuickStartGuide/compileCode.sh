#/bin/bash

# you can run it with command: bash compileCode.sh 1000 500 1e5 0

## compile 
cd build
# rm -rf CMakeCache.txt  CMakeFiles  cmake_install.cmake  Makefile
cmake ..
make
cd ..

if [[ -z "$1" ]]
then
    echo 'You did not define the number of particles'
    exit 0
fi

if [[ -z "$2" ]]
then
    echo 'You did not define the number of time steps'
    exit 0
fi

if [[ -z "$3" ]]
then
    echo 'You did not define delta t'
    exit 0
fi

if [[ -z "$4" ]]
then
    echo 'You did not define local diffusion'
    exit 0
fi

echo -e '\e[32m Number of particles you set is '$1'\e[0m'
echo -e '\e[33m Number of time steps you set is '$2'\e[0m'
echo -e '\e[35m Delta t you set is '$3'\e[0m'
echo -e '\e[36m Local diffusion you set is '$4'\e[0m'
echo '------------------------------------------'
echo '------------------------------------------'
echo '------------------------------------------'

./bin/main $1 $2 $3 $4
wait
./bin/Transform2DH5ParticleDataTo3D DFN_mesh.h5



