#!/bin/bash

# before push the code to github, run this bash script to clean the directory of the code

myArray=("BenchmarkCases" "ParticleTrackingFourFrac" "ClearBadH5FileForPercolationTest" \
  "CompareCPUGPU" "CompareWithAndWithoutDeadEnds" "ContinueParticleTransport" "FracturesFixedChangeDomainSize" \
  "ParticleTransportValidationSingleFrac" "PercolationTest" \
  "ReadMatAndOutputH5" "ReadObjectsAndDebug" "TestFunctionsOnlyForDevlopment" "TestResolutionEffect" \
  "Transform2DH5ParticleDataTo3D" "GenMultipleFamilies" "ConnectivityMultipleFamilies"
  "ColumnLikeDomainForDispersion" "ValidateLocalDiffusion" "ConvertH5Precision" "DispersionAtOneDensityValue" "IdenticalDFNDifferentPe" \
  "Dispersion_Rho_increasing" "ValidationOfDispersionPhDThesis")

for str in ${myArray[@]}; do
  cd ./SomeApplications/$str
  gio trash -f ./*.h5 ./*.m ./*.py
  gio trash -f ParticlePositionResult
  # cd build
  # gio trash -f CMakeCache.txt CMakeFiles cmake_install.cmake
  # cd ..
  # cd bin
  gio trash -f ./main* ./Transform2DH5ParticleDataTo3D ./ColumnLikeDomainForDispersion ./QuickStartGuide ./QuickStartGuide_* ./ContinueParticleTransport ./ValidateLocalDiffusion ./ConvertH5Precision
  cd ../../
  echo "---------------cleaned "$str
done

myArray_II=("QuickStartGuide" "CompilationTest" "bin")
for str in ${myArray_II[@]}; do
  cd  $str
  gio trash -f ./*.h5 ./*.m ./*.py
  gio trash -f ParticlePositionResult
  # cd build
  # gio trash -f CMakeCache.txt CMakeFiles cmake_install.cmake
  # cd ..
  # cd bin
  gio trash -f ./cuDFNsys_exe ./main* ./Transform2DH5ParticleDataTo3D ./ColumnLikeDomainForDispersion ./QuickStartGuide ./QuickStartGuide_* ./ContinueParticleTransport ./ValidateLocalDiffusion ./ConvertH5Precision
  cd ../
  echo "---------------cleaned "$str
done

cd ./lib
gio trash -f lib*.a
cd build 
gio trash -f CMakeCache.txt CMakeFiles cmake_install.cmake
cd ..
cd ..

echo "---------------cleaned lib"
