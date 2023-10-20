#/bin/bash

cd build 
cmake -DCMAKE_CUDA_ARCHITECTURES=60 ..
make
cd ..