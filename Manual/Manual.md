# Manual for _cuDFNsys_

* Name: Tingchang YIN

* Institution: Westlake University, China

* Date: Oct. 19, 2023

* Update date: Jan. 22, 2024 

* Email: yintingchang@foxmail.com

## 1. Compilation of a user interface code with _cuDFNsys_

Suppose that the _cuDFNsys_ library has been installed in `~/cuDFNsys/lib`, then `libcuDFNsys.a` can be seen in `~/cuDFNsys/lib`, which is a static library. How can one compile and install _cuDFNsys_ library (with nessesary dependencies, e.g., gmsh) can be seen in `README.md`.

Now, we can create a directory `CompilationTest` under `~/cuDFNsys`. We should test the compilation under this directory:
```
cd ~/cuDFNsys
mkdir CompilationTest
cd CompilationTest
touch Makefile
mkdir src
cd src
touch CompilationTest.cu
```

We can copy the following code to `~/cuDFNsys/CompilationTest/src/CompilationTest.cu`:
```
#include "cuDFNsys.cuh"
int main()
{
    int dev = 0;
    GPUErrCheck(cudaSetDevice(dev));
    cuDFNsys::Vector4<double> vec4 = cuDFNsys::MakeVector4(1.5, 1., 15., 0.);
    std::cout << "cuDFNsys::Vector4: " << vec4.x << ", " << vec4.y << ", " << vec4.z << ", " << vec4.w << ", " << std::endl;
    return 0;
}
```
The above code just generates a vector (type: `cuDFNsys::Vector4<double>`) named `vec4` of four elements (double precision), which can be accessed by `vec4.x`, `vec4.y`, `vec4.z` and `vec4.w`. This vector is actually `double3` in cuda!!!

To compile this code linking to `libcuDFNsys.a`, we can copy the following `Make` script to `~/cuDFNsys/CompilationTest/Makefile`:
```
ExeName=CompilationTest
# NVCC path
NVCC=/usr/lib/nvidia-cuda-toolkit/bin/nvcc
# include paths for headers
cuDFNsysIncludePath=$(HOME)/cuDFNsys/include
Hdf5IncludePath=/usr/include/hdf5/serial
GmshIncludePath=usr/include
EigenIncludePath=/usr/include
UmfpackIncludePath=/usr/include/suitesparse
# library paths
cuDFNsysLibraryPath=$(HOME)/cuDFNsys/lib
GmshLibraryPath=/usr/lib/x86_64-linux-gnu
UmfpackLibraryPath=/usr/lib/x86_64-linux-gnu
Hdf5LibraryPath=/usr/lib/x86_64-linux-gnu/hdf5/serial

INCDIRS=-I $(cuDFNsysIncludePath) \
		-I $(Hdf5IncludePath) \
		-I $(GmshIncludePath) \
		-I $(EigenIncludePath) \
		-I /usr/local/include \
		-I $(UmfpackIncludePath)
Lib_= -Xcompiler=-fopenmp \
	-lpthread \
	-lm \
	-lcudadevrt \
	-lcudart_static \
	-lrt \
	-lpthread \
	-ldl \
	-lgmsh \
	-lumfpack \
	-lamd \
	-lhdf5_cpp \
	-lhdf5 \
	-lsz \
	-lz \
	-ldl \
	-L$(GmshLibraryPath) \
	-L$(UmfpackLibraryPath) \
	-L$(Hdf5LibraryPath)
all: 
	$(NVCC) -DUSE_DOUBLES ./src/$(ExeName).cu $(cuDFNsysLibraryPath)/libcuDFNsys.a \
	-o ./main $(INCDIRS) $(Lib_) \
	-Xcudafe --display_error_number -Xcudafe --diag_suppress=3057 \
	-Xcudafe --diag_suppress=1301 \
	-Xcudafe --diag_suppress=3059 \
	-Xcudafe --diag_suppress=3060 \
	-Xcudafe --diag_suppress=3058 \
	-arch=sm_60 -std=c++17 -rdc=true
clean:
```

Note that, if necessary, some variables of this script should be changed for your Ubuntu, namely,
```
# NVCC path
NVCC=/usr/lib/nvidia-cuda-toolkit/bin/nvcc
# include paths for headers
cuDFNsysIncludePath=$(HOME)/cuDFNsys/include
Hdf5IncludePath=/usr/include/hdf5/serial
GmshIncludePath=usr/include
EigenIncludePath=/usr/include
UmfpackIncludePath=/usr/include/suitesparse
# library paths
cuDFNsysLibraryPath=$(HOME)/cuDFNsys/lib
GmshLibraryPath=/usr/lib/x86_64-linux-gnu
UmfpackLibraryPath=/usr/lib/x86_64-linux-gnu
Hdf5LibraryPath=/usr/lib/x86_64-linux-gnu/hdf5/serial
```
The above variables are the paths to different tools, CXX head files, and CXX libraries that _cuDFNsys_ depends. For example, you can change the variable `GmshIncludePath` like
```
GmshIncludePath=$(HOME)/gmsh/include
```
if you installed `gmsh-dev` under `~/gmsh` in your Ubuntu.

Now, we can comile `CompilationTest.cu` by:
```
cd ~/cuDFNsys/CompilationTest
make
```
If no error happens, then an executable file `main` should appear in `~/cuDFNsys/CompilationTest`. We just run it by
```
./main
```
and we should see
```
cuDFNsys::Vector4: 1.5, 1, 15, 0,
```
printed on the terminal window.

If errors happen when you run `./main`, e.g., `error while loading shared libraries: libgmsh.so: cannot open shared object file: No such file or directory`, two enviromental variables should be added to `~/.bashrc` as follows:
```
export LD_LIBRARY_PATH=path-to-gmsh-library:$LD_LIBRARY_PATH
export LIBRARY_PATH=path-to-gmsh-library:$LIBRARY_PATH
# change path-to-gmsh-library to the dynamic gmsh library in the computer
```
then update `~/.bashrc` by
```
source ~/.bashrc
```
This is a runtime error, that the computer does not know where the library `gmsh` is when running `./main`.

## 2. Use an executable file
To users who do not want to write a `.cpp` file, an executable file - `cuDFNsys_exe`, can be used.

First, the executable should be compiled.
```
cd ~/cuDFNsys/bin
make
```

Then, you can see two files, `cuDFNsys_exe` and `Transform2DH5ParticleDataTo3D`, generated in the current directory.

`cuDFNsys_exe` needs two input `.csv` files to get fracture attributes and other parameters.

To generate DFNs, I provide two `.csv` files.

### 2.1 csv to generate a deterministic DFN
```
#InputParametersForDeterministicDFN.csv
IfStochastic,0,
DomainSizeX,30,
DomainDimensionRatio,1,1,1,,,,,
Percolation_direction,2,,,,,,,
NumFractures,4,
Fracture_1, 1,0,0,0,0,0,10,0.2,5e-4,0,1,1,
Fracture_2, 0,0,1,0,0,0,55,0.1,7e-9,,,,
Fracture_3, 0,1,0,0,0,0,55,0.3,6e-7,
Fracture_4, 1,1,1,10,10,-25,4,0.3,6e-7,

```
This one means a deterministic DFN, defined by users.
`IfStochastic,0,` denotes that this DFN is deterministic.
`DomainSizeX,30,` denotes that the length of the side in the `x` direction is 30 m.
`DomainDimensionRatio,1,1,1,,,,,` denotes that the dimension ratio of the DFN domain is 1 : 1 : 1, thus the size of the domain is 30m, 30m, 30m
`Percolation_direction,2,,,,,,,` means that the program will check if the DFN is percolative in `z` direction
`NumFractures,4,` means that there are four fractures
`Fracture_1, 1,0,0,0,0,0,10,0.2,5e-4,0,1,1,`: there are 12 numbers separated by commas, representing the fracture attribues of a fracture. The first three `1,0,0` means the normal vector of the fracture plane. Then, the next three `0,0,0` is the center coodinate. The seventh `10` is the radius $R$ of the fracture. The eighth `0.2` is the $\beta$ value of fractures, the ninth `5e-4` is the $\gamma$ value. In _cuDFNsys_, the conductivity $k_f = b_f^3 / 12 =  \left[(\gamma R ^\beta) \right] ^ 3 / 12$, where $b_f$ is the aperture of fractures. The last three `0,1,1` denotes the coordinates of one vertices of the square fracture. If the last three number is `0,0,0,` or `,,,` (no numbers), then the program will randomly determine a vertex in the fracture plane with $R$.

### 2.2 csv to generate a random DFN
```
##InputParametersForStochasticDFN.csv
IfStochastic,1,
DomainSizeX,30,
DomainDimensionRatio,1,1,2,,,,,
Percolation_direction,2,,,,,,,
NumFractureGroups,3,,,,,,,
NumFractureEachGroup,70,50,30,,,,,
KappaValues,20,20,25,,,,,
MeanOrientationOfFisherDistribution,0,0,1,1,0,0,0,1,0,
ModeOfSizeDistribution,0,1,2,,,,,
SizeDistributionParameters,1.5,1,15,0,8.5,5.5,1,15,1,15,0,0,
Beta,0.2,0.3,0.4,
Gamma,2.0e-5,3.0e-6,4e-5,
```

`IfStochastic,1,` denotes a stochastic DFN.
`NumFractureGroups,3,,,,,,,` denotes that there three fracture groups.
`KappaValues,20,20,25,,,,,`: the $\kappa$ values (Fisher constants) of the three groups.
`MeanOrientationOfFisherDistribution,0,0,1,1,0,0,0,1,0,`: mean orientations of the three groups in vector forms.
`ModeOfSizeDistribution,0,1,2,,,,,`: the size distribution patterns of each fracture Group. When the mode is `0`, the size distribution is a power-law. When the mode  is `1`, the size distribution is lognormal. When the mode  is `2`, the size distribution is uniform. When the mode  is `3`, the size distribution is mono-sized.
`SizeDistributionParameters`: the statistic parameters of the size distributions of each group. Four numbers is a set of statistic parameters for a size distribution of a group. `1.5,1,15,0` means, as the distribution of first group is a power-law, the power-law exponent is 1.5, minimum $R$ is 1 m, maximum is 15 m. As the second group has a longnormal size distribution, `8.5,5.5,1,15` means that the mean is 8.5, variance is 5.5, min R is 1, max $R$ is 15 m. The this group has a uniform size distribution, thus `1,15,0,0,` means that min $R$ is 1, max $R$ is 15 m.
`Beta,0.2,0.3,0.4,` and `Gamma,2.0e-5,3.0e-6,4e-5,` mean the $\beta$ and $\gamma$ values of each group, respectively.

### 2.3. Run a DFN generation and mesh and flow simulation
One can run the program by `./cuDFNsys_exe`. Then the terminal shows
```
Creating a new DFN. Please provide the name of the .csv file, without the suffix .csv!
```
One can input the name of `.csv`, either `InputParametersForDeterministicDFN` or `InputParametersForStochasticDFN`.

Then, the program asks you to interactively input some other parameters, e.g., the minimum and maximum, the inlet and outlet head values, and so on.

### 2.4. Particle tracking
```
### PT_parameters.csv
NumParticles,20000,
NumTimeSteps,2000,
MolecularDiffusion,1.39715e-12,
DeltaT, 76491223.69119304,
InjectionMethod_initialcondition,Flux-weighted,
If_OutputAllPTInformationOrFPTCurve, 1,
SpacingOfControlPlanes, 30,
IfOutputVarianceOfDisplacementsEachStep, 1,
IfInjectAtCustomedPlane, 1,
CustomedPlaneInjection, 23,
IfUseFluxWeightedOrEqualProbableMixingIntersection, 1,
```
The parameters of particle tracking are defined in a `.csv`.

`NumParticles` is the number of particles injected into DFNs.
`InjectionMethod_initialcondition` is the initial condition, which can be `Flux-weighted`, `Resident`, or `Point`.
`If_OutputAllPTInformationOrFPTCurve` can be 0 or 1. The value `1` means particle coordinates of all steps are outputted.
`SpacingOfControlPlanes` is the spacing of control planes.
`IfOutputVarianceOfDisplacementsEachStep` can be 0 or 1. The value `1` means that the variance of particle displacements is outputted for each step.
`IfInjectAtCustomedPlane` can be 0 or 1. The value `1` means that the injection plane is user-defined.
`CustomedPlaneInjection`: if `IfInjectAtCustomedPlane == 1`, then `IfInjectAtCustomedPlane` is the location of injection plane.
`IfUseFluxWeightedOrEqualProbableMixingIntersection`: `0` is equiproble mixing rule at intersections. `1`: outgoing-flux-weighted mixing rule at intersections.

### 2.5. Visualization
After the simulation, one can see some `.m` or `.py` files in the current directory. One can run them by MATLAB or python to visualize DFN, DFN mesh and DFN hydraulic head field.

### 2.6. Data of DFN, mesh, and flow
These will be discussed in the next section.

## 3. Data files
By running `cuDFNsys_exe`, three `.h5` files are outputted, e.g., `Class_DFN.h5`, `Class_MESH.h5` and `Class_FLOW.h5`.

One can open them to see DFN parameters.

## 4. Quickstart examples: write a main.cpp

## 5. How to re-run a DFN


