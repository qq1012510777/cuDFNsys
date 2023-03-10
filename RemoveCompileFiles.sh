#!/bin/bash

# before push the code to github, run this bash script to clean the directory of the code

myArray=("BenchmarkCases" "ParticleTrackingFourFrac" "ClearBadH5FileForPercolationTest" \
  "CompareCPUGPU" "CompareWithAndWithoutDeadEnds" "ContinueParticleTransport" "FracturesFixedChangeDomainSize" \
  "ParticleTransportValidationSingleFrac" "PercolationTest" "QuickStartGuide" \
  "ReadMatAndOutputH5" "ReadObjectsAndDebug" "TestFunctionsOnlyForDevlopment" "TestResolutionEffect" \
  "Transform2DH5ParticleDataTo3D" "GenMultipleFamilies" "ConnectivityMultipleFamilies"
  "ColumnLikeDomainForDispersion")

for str in ${myArray[@]}; do
  cd $str
  gio trash -f ./*.h5 ./*.m ./*.py
  gio trash -f ParticlePositionResult
  cd build
  gio trash -f CMakeCache.txt CMakeFiles cmake_install.cmake
  cd ..
  cd bin
  gio trash -f ./main* ./Transform2DH5ParticleDataTo3D ./ColumnLikeDomainForDispersion
  cd ..
  cd ..
  echo "---------------cleaned "$str
done
