# Manual for _cuDFNsys_

* Name: Tingchang YIN

* Institution: Westlake University, China

* Date: Oct. 19, 2023

* Update date: Oct. 23, 2023 

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
The above code just generates a vector (type: `cuDFNsys::Vector4<double>`) named `vec4` of four elements (double precision), which can be accessed by `vec4.x`, `vec4.y`, `vec4.z` and `vec4.w`. More details about `cuDFNsys::Vector4<double>` are discussed later. This vector is actually `double3` in cuda.

To compile this code with `libcuDFNsys.a`, we can copy the following `Make` script to `~/cuDFNsys/CompilationTest/Makefile`:
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

Note that some variables of this script should be changed for your Ubuntu, namely,
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
if you installed `gmsh` under `~/gmsh` in your Ubuntu.

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

## 2. Documented _cuDFNsys_ classes 

_cuDFNsys_ is an object-oriented library. Flow and transport in DFNs can be simulated by **four** `CUDA C++` classes, i.e. `cuDFNsys::DFN<Data_type>`, `cuDFNsys::MeshDFN<Data_type>`, `cuDFNsys::FlowDFN<Data_type>` and `cuDFNsys::PTDFN<Data_type>`. The `Data_type` can be `double` or `float`, but the latter is not recommended, because it may lead to wrong results. The header file of the four classes is in `cuDFNsys/include/cuDFNsys.cuh`. The source file is in `cuDFNsys/src/cuDFNsys.cu`.

To use these classes, the steps are (1) establish an empty class; (2) setup the member variables; and (3) call member functions. 

### 2.1. DFN generation
To generate a DFN, some input parameters are required. For instance,two groups of fractures $g_1$ and $g_2$ are expected with
* Domain. A cuboid of $30 \times 30 \times 60$m.
* Fracture orientations. ($g_1$): $\kappa = 20$. ($g_2$): $\kappa = 10$. The orientations are described by the Fisher distribution with a Fisher constant $\kappa$. The mean orientations for the two groups are $(0, 0, 1)$ and $(1, 0, 0)$, respectively.
* Fracture sizes. Note that in _cuDFNsys_, the size of fractures is the radius $R$ of circumscribed circles of square fractures. ($g_1$): power-law, where the exponent $\alpha$ is 1.5, the minimum $R$ is 1, maximum $R$ is 15. ($g_2$): lognormal: the mean of logarithmic values is 8.5, the standard deviation of logarithmic values is 5.5, the minimum $R$ is 1, the maximum $R$ is 15.
* Fracture conductivity. In _cuDFNsys_, the conductivity $k_f$ of fractures is related to $R$ and the aperture $b$ of fractures. The relationship is $k_f = b^3 / 12 =  \left[(\gamma R ^\beta) \right] ^ 3 / 12$, where $\gamma$ and $\beta$ are a constant and an exponent, respectively. For $g_1$, $\gamma = 2\times 10^{-5}$, $\beta = 0.2$. For $g_2$, $\gamma = 3\times 10^{-5}$, $\beta = 0.3$.
* Number of fractures. ($g_1$): 70. ($g_2$): 80.

The _cuDFNsys_ class to generate a DFN is `cuDFNsys::DFN<T>` where `T` is a template (`double` or `float`).  The member variables includes
* `std::vector<int> NumFractures`: it is a `std::vector` of integers, denoting the number of fractures in each group
* `std::vector<T> Kappa`: it is a `std::vector` of `double`/`float` (depending on `T`), denoting the $\kappa$ value for each group.
* `std::vector<cuDFNsys::Vector3<T>> MeanOrientationOfFisherDistribution`: it is a `std::vector` of `double3`/`float3`, denoting the mean orientation for each group.
* `T DomainSizeX`: it is a `double`/`float`, denoteing the size of a domain in the `x` direction.
* `double3 DomainDimensionRatio`: it is a `double3`, denoteing the dimension ratio of the domain.
* `std::vector<T> Beta`: it is a `std::vector` of `double`/`float`, denoting the $\beta$ value for each group.
* `std::vector<T> Gamma`: it is a `std::vector` of `double`/`float`, denoting the $\gamma$ value for each group.
* `std::vector<int> ModeOfSizeDistribution`: it is a `std::vector` of integers, denoting the distribution pattern of fracture sizes for each group. When one element `= 0`, the size distribution is a power-law. When one element `= 1`, the size distribution is lognormal. When one element `= 2`, the size distribution is uniform. When one element `= 3`, the size distribution is mono-sized.
* `std::vector<cuDFNsys::Vector4<T>> SizeDistributionParameters`: it is a `std::vector` of `double4` or `float4`. For instance, there is only one fracture group. When `ModeOfSizeDistribution[0] = 0`, the size distribution is a power-law. `SizeDistributionParameters[0].x` is $\alpha$, `SizeDistributionParameters[0].y` is the minimum $R$, `SizeDistributionParameters[0].z` is the maximum $R$, `SizeDistributionParameters[0].w` means nothing. When `ModeOfSizeDistribution[0] = 1`, the size distribution is lognormal `SizeDistributionParameters[0].x` is the mean of logarithmic values, `SizeDistributionParameters[0].y` is the standard deviation of logarithmic values, `SizeDistributionParameters[0].z` is the minimum $R$, `SizeDistributionParameters[0].w` is the maximum $R$. When `ModeOfSizeDistribution[0] = 2`, the size distribution is uniform. `SizeDistributionParameters[0].x` is the minimum $R$, `SizeDistributionParameters[0].y` is the maximum $R$, `SizeDistributionParameters[0].z` and `SizeDistributionParameters[0].w` mean nothing. When `ModeOfSizeDistribution[0] = 3`, the size distribution is mono-sized. `SizeDistributionParameters[0].x` is $R$, `SizeDistributionParameters[0].y`, `SizeDistributionParameters[0].z` and `SizeDistributionParameters[0].w` mean nothing.

* `int PercoDir`: it is a integer in the range of $[0, 2]$, denoting the pre-defined percolation direction to be checked. 0: $x$ direction; 1: $y$ direction. 2: $z$ direction.
* `int NumFracturesTotal`: it is a integer, denoting the total number of fractures in the DFN.
* `unsigned long RandomSeed`: it is a unsigned long integer, denoting the random seed to be used to generated stochastic fractures.
* `thrust::host_vector<cuDFNsys::Fracture<T>> FracturesHost`: a `thrust::host_vector` of the struct `cuDFNsys::Fracture<T>`, denoting the data of fractures. The struct `cuDFNsys::Fracture<T>` will be described later.
* `thrust::device_vector<cuDFNsys::Fracture<T>> FracturesDevice`: a `thrust::device_vector` of the struct `cuDFNsys::Fracture<T>`, denoting the data of fractures. The struct `cuDFNsys::Fracture<T>` will be described later.
* `cuDFNsys::Fracture<T> *FracturesDevicePtr`: a pointer to `thrust::device_vector<cuDFNsys::Fracture<T>> FracturesDevice`
* `std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> IntersectionMap`: a `C++ `map. The key of `IntersectionMap` is a pair of `size_t`, which is the pair of IDs of fractures that intersect with each other. The value of the map is also a pair, which is a pair of `double3`/`float3`, denoting the coordinates of two ends of the intersections.
* `std::vector<std::vector<size_t>> ListClusters`: it is a `std::vector` of `std::vector`. It will store the clusters in DFNs. For instance, `ListClusters[0]`, `ListClusters[1]` and so on are clusters in a DFN. `ListClusters[0][0]`, `ListClusters[0][1]` and so on are IDs of fractures that belong to the same cluster.
* `std::vector<size_t> PercolationCluster`: a `std::vector` of `size_t`. It will store the ID of the CLUSTERs that are percolative. For example, if `Percolation_cluster[0] = 0`, it means `ListClusters[0]` is a percolative cluster.

The member functions of `cuDFNsys::DFN<T>`:
```
DFN(){};
void FractureGeneration();
void IdentifyIntersectionsClusters(const bool &IfTruncatedFractures);
void Visualization(const string &MatlabScriptName,
                   const string &PythonScriptName,
                   const string &HDF5FileName,
                   const bool &IfShowTruncatedFractures,
                   const bool &IfShowIntersections,
                   const bool &IfHightlighAllClusters,
                   const bool &IfShowOrientationDistribution);
void StoreInH5(const string &ClassNameH5);
void LoadClassFromH5(const string &ClassNameH5);
```
* `DFN` is an empty constructor of `cuDFNsys::DFN<T>`.
* `FractureGeneration` is called to generate fractures based on the member variables (if they are set properly). It does not require any inputs. This function will assign values to `FracturesHost` and `FracturesDevice`.
* `IdentifyIntersectionsClusters` is a function to be called to identify intersections of fractures. Intersections are stored in `IntersectionMap`.
* `Visualization` is called to do visualizations. Input parameters are required. `MatlabScriptName`: name of the matlab script, without the suffix `.m`. `PythonScriptName`: name of the python script. `HDF5FileName`: name of the HDF5 file storing the data. `IfShowTruncatedFractures`: if show the truncated fractures. `IfShowIntersections`: if show intersections. `IfHightlighAllClusters`: if highlight clusters in different color. `IfShowOrientationDistribution`: if show the distribution of orientations.
* `StoreInH5` is used to store the class data in a HDF5 file. `ClassNameH5`: the name of the HDF5 file.
* `LoadClassFromH5`: is used to generate a class by loading data from a HDF5 file (generated by `StoreInH5`).

Examples of `cuDFNsys::DFN<T>` can be seen in `cuDFNsys/QuickStartGuide/src/QuickStartGuide.cu` and `cuDFNsys/QuickStartGuide/src/QuickStartGuide_DFN_I_DFN.cu`.

In `cuDFNsys/QuickStartGuide/src/QuickStartGuide.cu`, the relevant code is:
```
time_t t;
time(&t);
cuDFNsys::DFN<double> my_dfn;
my_dfn.NumFractures = {70, 80};
my_dfn.Kappa = {20, 10};
my_dfn.MeanOrientationOfFisherDistribution = {make_double3(0., 0., 1.),
                                              make_double3(1., 0., 0.)};
my_dfn.DomainSizeX = 30;
my_dfn.DomainDimensionRatio = make_double3(1., 1., 2.);
my_dfn.Beta = {0.2, 0.3};
my_dfn.Gamma = {2.0e-5, 3.0e-6};
my_dfn.ModeOfSizeDistribution = {0, 1};
my_dfn.SizeDistributionParameters = {make_double4(1.5, 1., 15., 0.),
                                     make_double4(8.5, 5.5, 1., 15.)};
my_dfn.PercoDir = 2;
my_dfn.RandomSeed = (unsigned long)t;
my_dfn.FractureGeneration();
my_dfn.IdentifyIntersectionsClusters(true);
my_dfn.Visualization("DFN_VISUAL", "DFN_VISUAL", "DFN_VISUAL", true, true,
                     true, true);
```
In this code, `t` is generated as a random seed. An empty class of `cuDFNsys::DFN<double>` is firstly generated. Then the member variables, e.g., `NumFractures` and `RandomSeed`, are set. Those variables match with the parameters stated in the first paragraph in this subsection. Then, member functions `FractureGeneration`, `IdentifyIntersectionsClusters` and `Visualization` are called in order. Visualization scripts, e.g., `DFN_VISUAL.m` and `DFN_VISUAL.py` are generated in the current working directory. We can see the DFN by running the `DFN_VISUAL.m` is MATLAB, or by running the command `python3 DFN_VISUAL.py` in Ubuntu terminal (`mayavi` is required). One DFN is shown below.

<p align="center">
  <img width="500" src="https://github.com/qq1012510777/cuDFNsys/blob/main/Manual/A_DFN.png">
</p>
<p align="center">
    <em>Visualization of a DFN. The color and color bar denote the $z$ value of fractures. </em>
</p>

In `cuDFNsys/QuickStartGuide/src/QuickStartGuide_DFN_I_DFN.cu`, the code is similar, but it stores the class data in a HDF5 file:
```
time_t t;
time(&t);
cuDFNsys::DFN<double> my_dfn;
my_dfn.NumFractures = {70, 80};
my_dfn.Kappa = {20, 10};
my_dfn.MeanOrientationOfFisherDistribution = {make_double3(0., 0., 1.),
                                              make_double3(1., 0., 0.)};
my_dfn.DomainSizeX = 30;
my_dfn.DomainDimensionRatio = make_double3(1., 1., 2.);
my_dfn.Beta = {0.2, 0.3};
my_dfn.Gamma = {2.0e-5, 3.0e-6};
my_dfn.ModeOfSizeDistribution = {0, 1};
my_dfn.SizeDistributionParameters = {make_double4(1.5, 1., 15., 0.),
                                     make_double4(8.5, 5.5, 1., 15.)};
my_dfn.PercoDir = 2;
my_dfn.RandomSeed = (unsigned long)t;
my_dfn.FractureGeneration();
my_dfn.IdentifyIntersectionsClusters(true);
my_dfn.Visualization("DFN_VISUAL_I", "DFN_VISUAL_I", "DFN_VISUAL_I", true,
                     true, true, true);
my_dfn.StoreInH5("Class_DFN");
```
`my_dfn.StoreInH5("Class_DFN")` outputs a HDF5 file named `Class_DFN.h5`. In `QuickStartGuide_DFN_II_MESH.cu`, the HDF5 file is loaded and then the data is assigned to a class of `cuDFNsys::DFN<double>`.

