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

_cuDFNsys_ is an object-oriented library. Flow and transport in DFNs can be simulated by just **four** `CUDA C++` classes, i.e. `cuDFNsys::DFN<Data_type>`, `cuDFNsys::MeshDFN<Data_type>`, `cuDFNsys::FlowDFN<Data_type>` and `cuDFNsys::PTDFN<Data_type>`. The `Data_type` can be `double` or `float`, but the latter is not recommended, because it may lead to wrong results. The header file of the four classes is in `cuDFNsys/include/cuDFNsys.cuh`. The source file is in `cuDFNsys/src/cuDFNsys.cu`.

To use these classes, the steps are (1) establish an empty class; (2) setup the member variables; and (3) call member functions. 

### 2.1. DFN generation
Fractrues in _cuDFNsys_ are 3D squares. The distribution of fracture centers is uniform right now.

To generate a DFN, some input parameters are required. For instance,two groups of fractures $g_1$ and $g_2$ are expected with
* **Domain**. A cuboid of $30 \times 30 \times 60$m. Hence the dimension ratio is $(1, 1, 2)$.
* **Fracture orientations**. ($g_1$): $\kappa = 20$. ($g_2$): $\kappa = 10$. The orientations are described by the Fisher distribution with a Fisher constant $\kappa$. The mean orientations for the two groups are $(0, 0, 1)$ and $(1, 0, 0)$, respectively. The mean orientation of the Fisher distribution means that the up-pointing unit normal vector of planar fractures will cluster about the mean orientation, and the clustering degree is denoted by $\kappa$.
* **Fracture sizes**. Note that in _cuDFNsys_, the size of fractures is the radius $R$ of circumscribed circles of square fractures. ($g_1$): power-law, where the exponent $\alpha$ is 1.5, the minimum $R$ is 1, maximum $R$ is 15. ($g_2$): lognormal: the mean of logarithmic values is 8.5, the standard deviation of logarithmic values is 5.5, the minimum $R$ is 1, the maximum $R$ is 15.
* **Fracture conductivity**. In _cuDFNsys_, the conductivity $k_f$ of fractures is related to $R$ and the aperture $b$ of fractures. The relationship is $k_f = b^3 / 12 =  \left[(\gamma R ^\beta) \right] ^ 3 / 12$, where $\gamma$ and $\beta$ are a constant and an exponent, respectively. For $g_1$, $\gamma = 2\times 10^{-5}$, $\beta = 0.2$. For $g_2$, $\gamma = 3\times 10^{-6}$, $\beta = 0.3$.
* **Number of fractures**. ($g_1$): 70. ($g_2$): 80.

#### 2.1.1 Member variables of `cuDFNsys::DFN`

The _cuDFNsys_ class to generate a DFN is `cuDFNsys::DFN<T>` where `T` is a template (`double` or `float`).  The member variables includes
* **`std::vector<int> NumFractures`**: it is a `std::vector` of integers, denoting the number of fractures in each group. So the number of fracture groups is `NumFractures.size()`.
* **`std::vector<T> Kappa`**: it is a `std::vector` of `double`/`float` (depending on `T`), denoting the $\kappa$ value for each group.
* **`std::vector<cuDFNsys::Vector3<T>> MeanOrientationOfFisherDistribution`**: it is a `std::vector` of `double3`/`float3`, denoting the mean orientation for each group.
* **`T DomainSizeX`**: it is a `double`/`float`, denoting the size of a domain in the `x` direction.
* **`double3 DomainDimensionRatio`**: it is a `double3`, denoting the dimension ratio of the domain.
* **`std::vector<T> Beta`:** it is a `std::vector` of `double`/`float`, denoting the $\beta$ value for each group.
* **`std::vector<T> Gamma`**: it is a `std::vector` of `double`/`float`, denoting the $\gamma$ value for each group.
* **`std::vector<int> ModeOfSizeDistribution`**: it is a `std::vector` of integers, denoting the distribution pattern of fracture sizes for each group. When one element `= 0`, the size distribution is a power-law. When one element `= 1`, the size distribution is lognormal. When one element `= 2`, the size distribution is uniform. When one element `= 3`, the size distribution is mono-sized.
* **`std::vector<cuDFNsys::Vector4<T>> SizeDistributionParameters`**: it is a `std::vector` of `double4` or `float4`. For instance, there is only one fracture group. When `ModeOfSizeDistribution[0] = 0`, the size distribution is a power-law. `SizeDistributionParameters[0].x` is $\alpha$, `SizeDistributionParameters[0].y` is the minimum $R$, `SizeDistributionParameters[0].z` is the maximum $R$, `SizeDistributionParameters[0].w` means nothing. When `ModeOfSizeDistribution[0] = 1`, the size distribution is lognormal `SizeDistributionParameters[0].x` is the mean of logarithmic values, `SizeDistributionParameters[0].y` is the standard deviation of logarithmic values, `SizeDistributionParameters[0].z` is the minimum $R$, `SizeDistributionParameters[0].w` is the maximum $R$. When `ModeOfSizeDistribution[0] = 2`, the size distribution is uniform. `SizeDistributionParameters[0].x` is the minimum $R$, `SizeDistributionParameters[0].y` is the maximum $R$, `SizeDistributionParameters[0].z` and `SizeDistributionParameters[0].w` mean nothing. When `ModeOfSizeDistribution[0] = 3`, the size distribution is mono-sized. `SizeDistributionParameters[0].x` is $R$, `SizeDistributionParameters[0].y`, `SizeDistributionParameters[0].z` and `SizeDistributionParameters[0].w` mean nothing.
* **`int PercoDir`**: it is a integer in the range of $[0, 2]$, denoting the pre-defined percolation direction to be checked. 0: $x$ direction; 1: $y$ direction. 2: $z$ direction.
* **`unsigned long RandomSeed`**: it is a unsigned long integer, denoting the random seed to be used to generated stochastic fractures.
>[!NOTE]
>
> **The member variables listed above should be defined before the exact generation of a DFN**.
> **Fracture data and other information will be populated into the following member variables after some member functions are called**

* **`int NumFracturesTotal`**: it is a integer, denoting the total number of fractures in the DFN.
* **`thrust::host_vector<cuDFNsys::Fracture<T>> FracturesHost`**: a `thrust::host_vector` of the struct `cuDFNsys::Fracture<T>`, denoting the data of fractures. The struct `cuDFNsys::Fracture<T>` will be described later.
* **`thrust::device_vector<cuDFNsys::Fracture<T>> FracturesDevice`**: a `thrust::device_vector` of the struct `cuDFNsys::Fracture<T>`, denoting the data of fractures. The struct `cuDFNsys::Fracture<T>` will be described later.
* **`cuDFNsys::Fracture<T> *FracturesDevicePtr`**: a pointer to `thrust::device_vector<cuDFNsys::Fracture<T>> FracturesDevice`
* **`std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> IntersectionMap`**: a `C++ `map. The key of `IntersectionMap` is a pair of `size_t`, which is the pair of IDs of fractures that intersect with each other. The value of the map is also a pair, which is a pair of `double3`/`float3`, denoting the coordinates of two ends of the intersections.
* **`std::vector<std::vector<size_t>> ListClusters`**: it is a `std::vector` of `std::vector`. It will store the clusters in DFNs. For instance, `ListClusters[0]`, `ListClusters[1]` and so on are clusters in a DFN. `ListClusters[0][0]`, `ListClusters[0][1]` and so on are IDs of fractures that belong to the same cluster.
* **`std::vector<size_t> PercolationCluster`**: a `std::vector` of `size_t`. It will store the ID of the CLUSTERs that are percolative. For example, if `Percolation_cluster[0] = 0`, it means `ListClusters[0]` is a percolative cluster.

#### 2.1.2. Struct of a single fracture: `cuDFNsys::Fracture<T>`
`cuDFNsys::Fracture<T>` is a struct to represent a fracture. The member variables of this struct are
```
template <typename T>
struct Fracture
{
public:
    // conductivity of the fracture
    T Conductivity;

    // 3D verts
    cuDFNsys::Vector3<T> Verts3D[4];

    // center
    cuDFNsys::Vector3<T> Center;

    // verts of 3D square/fracture
    cuDFNsys::Vector3<T> Verts3DTruncated[8];

    // number of verts of truncated fractures
    int NumVertsTruncated;

    // radius of circumscribe circle
    T Radius;

    // if the fracture intersects the six faces of cubic model?
    bool ConnectModelSurf[6]; // x_min, x_max, y_min, y_max, z_min, z_max,

    // normal vec (normalized)
    cuDFNsys::Vector3<T> NormalVec;
}
```  
* **`Conductivity`** is $k_f$.
* **`Verts3D`** is the four vertices of the square fractures. 
* **`center`** is the center of the fracture. 
* **`Verts3DTruncated`** is the vertices of square fractures after the fracture is truncated by the domain. 
* **`NumVertsTruncated`** is the number of vertices after truncation, so its maximum number is eight.  
* **`Radius`** is the radius of fractures. 
* **`ConnectModelSurf`** denotes if the fracture intersects the six surfaces of the domain, `true` means intersected. 
* **`NormalVec`** is a `double3` or `float3`, denoting the normal vector (always pointing up) of a fracture. 

#### 2.1.3. Member functions of `cuDFNsys::DFN`
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
* **`DFN()`** is an empty constructor of `cuDFNsys::DFN<T>`.
* **`FractureGeneration()`** is called to generate fractures based on the member variables (if they are defined properly). It does not require any inputs. This function will generate a DFN and fracture data are populated into the vectors, i.e., `FracturesHost` and `FracturesDevice`.
* **`IdentifyIntersectionsClusters`** is a function to be called to identify intersections of fractures. Intersection data are stored in `IntersectionMap`.
* **`Visualization()`** is called to do visualizations. Input parameters are required. `MatlabScriptName`: name of the matlab script, without the suffix `.m`. `PythonScriptName`: name of the python script. `HDF5FileName`: name of the HDF5 file storing the data. `IfShowTruncatedFractures`: if show the truncated fractures. `IfShowIntersections`: if show intersections. `IfHightlighAllClusters`: if highlight clusters in different color. `IfShowOrientationDistribution`: if show the distribution of orientations.
* **`StoreInH5()`** is used to store the class data in a HDF5 file. `ClassNameH5`: the name of the HDF5 file.
* **`LoadClassFromH5()`**: is used to generate a class by loading data from a HDF5 file (generated by `StoreInH5`).

Examples of `cuDFNsys::DFN<T>` can be seen in `cuDFNsys/QuickStartGuide/src/QuickStartGuide.cu` and `cuDFNsys/QuickStartGuide/src/QuickStartGuide_DFN_I_DFN.cu`.

#### 2.1.4. Examples of `cuDFNsys::DFN`
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
In this code, `t` is generated as a random seed. An empty class of `cuDFNsys::DFN<double>` is firstly generated. Then the member variables, e.g., `NumFractures` and `RandomSeed`, are set. Those variables match with the parameters stated in the first paragraph in Section 2.1.1. Then, member functions `FractureGeneration`, `IdentifyIntersectionsClusters` and `Visualization` are called in order. Visualization scripts, e.g., `DFN_VISUAL.m` and `DFN_VISUAL.py` are generated in the current working directory. We can see the DFN by running the `DFN_VISUAL.m` is MATLAB, or by running the command `python3 DFN_VISUAL.py` in Ubuntu terminal (`mayavi` is required). One DFN is shown below.

<p align="center">
  <img width="500" src="https://github.com/qq1012510777/cuDFNsys/blob/main/Manual/A_DFN.png">
</p>
<p align="center">
    <em>Visualization of a DFN (by python-mayavi). The color and color bar denote the $z$ value of fractures. </em>
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
`my_dfn.StoreInH5("Class_DFN")` outputs a HDF5 file named `Class_DFN.h5`. In a following code `QuickStartGuide_DFN_II_MESH.cu`, you can see there that the HDF5 file is loaded and then the data are populated into an empty class of `cuDFNsys::DFN<double>`, then an identical DFN is generated without using of `FractureGeneration()`.

### 2.2. Mesh generation
In `cuDFNsys/QuickStartGuide/src/QuickStartGuide.cu`, the mesh generation is very simple, just several lines:
```
    cuDFNsys::MeshDFN<double> meshGen;
    meshGen.MinElementSize = 1;
    meshGen.MaxElementSize = 3;
    meshGen.MeshGeneration(my_dfn);
    meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL", "DFN_MESH_VISUAL",
                          "DFN_MESH_VISUAL", true, true);
```
Before running `meshGen.MeshGeneration()`, only the constrain for the size of finite element is set. Note that the GMSH C++ API is used to generate a conforming mesh in a DFN, and the constrain for the elements is just roughly satisfied, meaning that there would be a little bigger or smaller elements.

#### 2.2.1. Member variables of `cuDFNsys::MeshDFN`
Member variables of this class includes:
* **`T MinElementSize`**: it is a `double`/`float`, denoting the minimum element size you want.
* **`T MaxElementSize`**: it is a `double`/`float`, denoting the maximum element size you want.
* **`T MeanGridSize`**: it is a `double`/`float`, denoting the average size of elements after mesh.
* **`std::vector<size_t> FracsPercol`**: it is a vector of `size_t`, denoting the IDs of fractures belonging to the percolation cluster
* **`cuDFNsys::Mesh<T> MeshData`**: this is also a class, even it is a member of `cuDFNsys::MeshDFN`. It has some member variables to store the mesh data. It will be described later.

#### 2.2.2. Member functions of `cuDFNsys::MeshDFN`
The member functions of this class includes:
```
MeshDFN(){};
void MeshGeneration(cuDFNsys::DFN<T> &my_dfn);
void Visualization(cuDFNsys::DFN<T> my_dfn,
                   const string &MatlabScriptName,
                   const string &PythonScriptName,
                   const string &HDF5FileName,
                   const bool &IfCheck2DCoordinatesOfMesh,
                   const bool &IfCheckEdgeAttributes);
void StoreInH5(const string &ClassNameH5);
void LoadClassFromH5(const string &ClassNameH5);
```
* **`MeshDFN()`** is used to establish an empty class of `cuDFNsys::MeshDFN`.
* **`MeshGeneration()`** is used to generate a mesh. It requires `cuDFNsys::DFN` as inputs.
* **`Visualization()`** is called to do visualizations. Input parameters are required. `my_dfn` is a class of `cuDFNsys::DFN`. `MatlabScriptName`: name of the matlab script, without the suffix `.m`. `PythonScriptName`: name of the python script. `HDF5FileName`: name of the HDF5 file storing the data. `IfCheck2DCoordinatesOfMesh`: not very impotent, just set it to be `false`. `IfCheckEdgeAttributes`: highligh the mesh edges in different colors according to their attributes (at interior or boundary?)
* **`StoreInH5()`** is used to store the class data in a HDF5 file. `ClassNameH5`: the name of the HDF5 file.
* **`LoadClassFromH5()`**: is used to generate a class by loading data from a HDF5 file (generated by `StoreInH5`).

#### 2.2.3. Class of `cuDFNsys::Mesh<T>`
`cuDFNsys::Mesh<T> MeshData` is a member variable in `cuDFNsys::MeshDFN`, which is also a class. Here I show some important memeber variables of `cuDFNsys::Mesh<T>`:
```
template <typename T>
class Mesh
{
public:
    thrust::host_vector<cuDFNsys::Vector3<T>> Coordinate3D;
    thrust::host_vector<uint3> Element3D;
    thrust::host_vector<uint> ElementFracTag; // from 0
}
```
* **`Coordinate3D`**: it is a `thrust::host_vector` of `double3` or `float3`. Each element is the 3D coordinate of one node in the mesh.
* **`Element3D`**: it is a `thrust::host_vector` of `uint3`. Each element contains the three IDs of the nodes that consist of a finite element.
* **`ElementFracTag`**: it is a `thrust::host_vector` of `uint`. Ech element is the ID of fracture that the element lies on.

In this class, the GMSH C++ API is called to generate a mesh in a DFN, and the mesh is conforming.

After meshing, mesh data will be populated into these member variable. Therefore, by assessing these variables, the mesh information can be obtained, and they are used in the subsequent FEM analysis.

#### 2.2.4 Examples of `cuDFNsys::MeshDFN`
In `cuDFNsys/QuickStartGuide/src/QuickStartGuide.cu`, the relevant code is:
```
cuDFNsys::MeshDFN<double> meshGen;
meshGen.MinElementSize = 1;
meshGen.MaxElementSize = 3;
meshGen.MeshGeneration(my_dfn);
meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL", "DFN_MESH_VISUAL",
                      "DFN_MESH_VISUAL", true, true);
```
First, an empty class of `cuDFNsys::MeshDFN<double>` is generated and named `meshGen`. Then the constrain for element sizes is set by defining `meshGen.MinElementSize` and 
`meshGen.MaxElementSize`. Next, the mesh is generated by calling `meshGen.MeshGeneration()`. Finally, the visualization script and data are output by `meshGen.Visualization`. A DFN mesh is shown below.

<p align="center">
  <img width="500" src="https://github.com/qq1012510777/cuDFNsys/blob/main/Manual/A_DFN_MESH.png">
</p>
<p align="center">
    <em>Visualization of a DFN mesh (by python-mayavi). The color and color bar denote the $z$ value of finite elements. </em>
</p>

In `cuDFNsys/QuickStartGuide/src/QuickStartGuide_DFN_II_MESH.cu`, it loads a DFN from a HDF5 file, then generates a mesh, finally stores the mesh data in a HDF5 file:
```
cuDFNsys::DFN<double> my_dfn;
my_dfn.LoadClassFromH5("Class_DFN");
my_dfn.Visualization("DFN_VISUAL_II", "DFN_VISUAL_II", "DFN_VISUAL_II",
                     true, true, true, true);
cuDFNsys::MeshDFN<double> meshGen;
meshGen.MinElementSize = 1;
meshGen.MaxElementSize = 3;
meshGen.MeshGeneration(my_dfn);
meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL_I", "DFN_MESH_VISUAL_I",
                      "DFN_MESH_VISUAL_I", true, true);
meshGen.StoreInH5("Class_MESH");
my_dfn.StoreInH5("Class_DFN");
```
The DFN data is loaded by `my_dfn.LoadClassFromH5`. The mesh data is stored by `meshGen.StoreInH5`. Note that in the end, the DFN data is stored again, because in the mesh step, those fractures that do not contribute to the flow are removed from the `thrust::host_vector`. That is to say, the class `my_dfn` is actually updated after meshing.

### 2.3. Flow simulation
The flow is numerically solved by a mixed hybrid finite element method, in which the continuity of flux can be satisfied. Note that regular node-based finite element methods cannot achieve that. The flow simulation is implement by a class named `cuDFNsys::FlowDFN<T>` where `T` is a template (`double` is recommanded).

#### 2.3.1. Member variables of `cuDFNsys::FlowDFN`
The member variables of `cuDFNsys::FlowDFN` includes:
* **`T InletHead`**: it is a `double` or `float`, denoting the hydraulic head at the inlet plane.
* **`T OutletHead`**: it is a `double` or `float`, denoting the hydraulic head at the outlet plane.

>[!NOTE]
>
> **The member variables listed above should be defined before the flow simulation in a DFN**.
> **Flow velocity and other information will be populated into the following member variables after some member functions are called**

* **`T MaxVelocity`**: it is a `double` or `float`, denoting the maximum fracture velocity $[LT^{-1}]$ in a DFN.
* **`T MeanVelocity`**: it is a `double` or `float`, denoting the minimum fracture velocity $[LT^{-1}]$ in a DFN.
* **`cuDFNsys::MHFEM<T> FlowData`**: it is a class with member variables storing the flow data, e.g., velocity at each edge of elements. This class will be described later.

#### 2.3.2. Member functions of `cuDFNsys::FlowDFN`
The member functions includes:
```
FlowDFN(){};
void FlowSimulation(cuDFNsys::DFN<T> my_dfn,
                    cuDFNsys::MeshDFN<T> my_mesh);
void Visualization(cuDFNsys::DFN<T> my_dfn,
                   cuDFNsys::MeshDFN<T> my_mesh,
                   const string &MatlabScriptName,
                   const string &PythonScriptName,
                   const string &HDF5FileName);
void StoreInH5(const string &ClassNameH5);
void LoadClassFromH5(const string &ClassNameH5);
```
* **`FlowSimulation()`** is called to generated mesh, after `InletHead` and `OutletHead` are defined. It requires a `cuDFNsys::DFN<T>` and a `cuDFNsys::MeshDFN<T>` as inputs. 
* **`Visualization()`** is called to do visualizations. Input parameters are required. `my_dfn` is a class of `cuDFNsys::DFN`. `my_mesh` is a class of `cuDFNsys::MeshDFN`. `MatlabScriptName`: name of the matlab script, without the suffix `.m`. `PythonScriptName`: name of the python script. `HDF5FileName`: name of the HDF5 file storing the data.
* **`StoreInH5()`** is used to store the class data in a HDF5 file. `ClassNameH5`: the name of the HDF5 file.
* **`LoadClassFromH5()`**: is used to generate a class by loading data from a HDF5 file (generated by `StoreInH5`).

#### 2.3.3. Class of `cuDFNsys::MHFEM<T>`
`cuDFNsys::MHFEM<T>` is a class, which has some impotant variables:
```
template <typename T>
class MHFEM
{
public:
    Eigen::MatrixXd PressureInteriorEdge;
    Eigen::MatrixXd PressureEles;
    Eigen::MatrixXd VelocityNormalScalarSepEdges;
};
```
* **`Eigen::MatrixXd PressureInteriorEdge`**: it is a Eigen vector of `double`. Each element is the hydraulic head value of an interior edge.
* **`Eigen::MatrixXd PressureEles`**: it is a Eigen vector of `double`. Each element is the average hydraulic head value of a finite element of the mesh.
* **`Eigen::MatrixXd VelocityNormalScalarSepEdges`**: it is a Eigen vector of `double`. Each element is the normal volumetric flux rate per width $\left[L^2T^{-1}\right]$ of the edge of a finite element.

#### 2.3.4. Examples of `cuDFNsys::FlowDFN`
In `cuDFNsys/QuickStartGuide/src/QuickStartGuide.cu`, the relevant code is:
```
cuDFNsys::FlowDFN<double> flowDFN;
flowDFN.InletHead = 60;
flowDFN.OutletHead = 0;
flowDFN.FlowSimulation(my_dfn, meshGen);
flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL", "DFN_FLOW_VISUAL",
                          "DFN_FLOW_VISUAL");
``` 
Firstly, an empty class of `cuDFNsys::FlowDFN<double>` is created. Then the inlet and outlet hydralic head values are set. Next, the flow is solved by calling `flowDFN.FlowSimulation()`. Finally, visualization scripts are output by `flowDFN.Visualization()`. An example of flow simulation is shown below.

<p align="center">
  <img width="500" src="https://github.com/qq1012510777/cuDFNsys/blob/main/Manual/A_DFN_FLOW.png">
</p>
<p align="center">
    <em>Visualization of the flow in a DFN (by python-mayavi). The color and color bar denote the head values. </em>
</p>

In `cuDFNsys/QuickStartGuide/src/QuickStartGuide_DFN_III_FLOW.cu`, it loads a DFN and the mesh data from HDF5 files, then solves the steady-state flow, finally stores the flow data in a HDF5 file:
```
cuDFNsys::DFN<double> my_dfn;
my_dfn.LoadClassFromH5("Class_DFN");
cuDFNsys::MeshDFN<double> meshGen;
meshGen.LoadClassFromH5("Class_MESH");
meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL_II", "DFN_MESH_VISUAL_II",
                      "DFN_MESH_VISUAL_II", true, true);
cuDFNsys::FlowDFN<double> flowDFN;
flowDFN.InletHead = 60;
flowDFN.OutletHead = 0;
flowDFN.FlowSimulation(my_dfn, meshGen);
flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL_I",
                      "DFN_FLOW_VISUAL_I", "DFN_FLOW_VISUAL_I");
flowDFN.StoreInH5("Class_FLOW");
```
The DFN data is loaded by `my_dfn.LoadClassFromH5`. The mesh data is loaded by `meshGen.LoadClassFromH5`.

### 2.4. Particle tracking
The class `cuDFNsys::PTDFN` is used to run particle tracking in a DFN. This class requires a lot of set-ups. The Peclet number $Pe$ is defined as
$$Pe = \frac{V_c \times l_c}{D_m}$$
where $V_c$ is the characteristic length scale in DFNs, $l_c$ is the characteristic length scale, and $D_m$ is the molecular diffusion coefficient for random walkers. Before the particle tracking, the values of $Pe$, $V_c$ and $l_c$ should be pre-defined, and $D_m$ can be then determined. In the code `QuickStartGuide.cu` or `QuickStartGuide_DFN_IV_PT.cu`, $V_c$ is set to be the average fracture velocity.

#### 2.4.1. Member variables of `cuDFNsys::PTDFN`
The member variables of `cuDFNsys::PTDFN` includes:
* **`int NumParticles`**: it is the number of particles that will be injected into the DFN.
* **`int NumTimeSteps`**: it is the number of time steps one wants to run.
* **`T PecletNumber`**: it is a `double`/`float`, denoting the value of $Pe$.
* **`T LengthScalePe`**: it is a `double`/`float`, denoting the value of $l_c$.
* **`T VelocityScalePe`**: it is a `double`/`float`, denoting the value of $V_c$.
* **`T MolecularDiffusion`**: it is a `double`/`float`, denoting the value of $D_m$.
* **`T TimeScaleCrossElement`**: it is a `double`/`float`, denoting that the time scale in a DFN for a random walker to cross a finite element.
* **`T FactorTimeScaleCrossElement`**: it is a `double`/`float`, which is a reduction factor to reduce the `TimeScaleCrossElement`.
* **`T DeltaT`**: it is a `double`/`float`, denoting $\delta t$ for a time step.
* **`bool FluexWeightedOrUniformInjection`**: it is a `bool`. `true` means the use of flux-weighted injection. `false` means the use of resident injection.
* **`bool IfUseFluxWeightedOrEqualProbableMixingIntersection`**: it is a `bool`. `true` means the use of outgoing-flux-weighted mixing rule at intersections. `false` means the use of equiprobable mixing rule at intersections.
* **`T SpacingOfControlPlanes`**: it is a `double` or `float`, denoting the spacing of control planes in the DFN. A control plane means that once random particles is cross the plane for the first time, the travel time and other information will be recorded.
* **`bool IfOutputVarianceOfDisplacementsEachStep`**: it is a `bool`. `true` means that the variance of solute particle displacements are recorded for each time step.
* **`bool IfInjectAtCustomedPlane`**: it is a `bool`. `true` means that the injection plane is not at the inlet, but at `CustomedPlaneInjection`.
* **`T CustomedPlaneInjection`**: it is a `double` or `float`. By using it, one can set the injection plane at any location.
* **`bool OutputAllPTInformationOrFPTCurve`**: it is a `bool`. `true` means that the particles' coodinates for each step will be stored in HDF5 files. `false` means that only the spatial information at last step is output, and the first passage times are recorded.

#### 2.4.2. Member functions of `cuDFNsys::PTDFN`
The member functions of `cuDFNsys::PTDFN` includes:
```
PTDFN(){};
void ParticleTracking(cuDFNsys::DFN<T> my_dfn,
                      cuDFNsys::MeshDFN<T> my_mesh,
                      cuDFNsys::FlowDFN<T> my_flow);
void Visualization(cuDFNsys::DFN<T> my_dfn,
                   cuDFNsys::MeshDFN<T> my_mesh,
                   cuDFNsys::FlowDFN<T> my_flow,
                   const string &MatlabScriptName,
                   const string &PythonScriptName,
                   const string &HDF5FileNameOfFlowDFN);
```
* **`PTDFN()`**: it is an empty constructor to generate an empty class of `cuDFNsys::PTDFN`
* **`ParticleTracking()`** is called to do particle tracking. It required some inputs: `my_dfn` is a class of `cuDFNsys::DFN`, `my_mesh` is a class of `cuDFNsys::MeshDFN<T>`, and `my_flow` is a class of `cuDFNsys::FlowDFN<T> `.
* **`Visualization()`** is called to do visualizations. Input parameters are required. `my_dfn` is a class of `cuDFNsys::DFN`. `my_mesh` is a class of `cuDFNsys::MeshDFN`, `my_flow` is a class of `cuDFNsys::FlowDFN<T>`. `MatlabScriptName`: name of the matlab script, without the suffix `.m`. `PythonScriptName`: name of the python script. `HDF5FileName`: name of the HDF5 file storing the visualization of **FLOW** data.

#### 2.4.3. Examples of `cuDFNsys::PTDFN`
In the code `cuDFNsys/QuickStartGuide/src/QuickStartGuide.cu`, the relevant code is:
```
cuDFNsys::PTDFN<double> particleTracking;
particleTracking.NumParticles = 20000;
particleTracking.NumTimeSteps = 200;
particleTracking.PecletNumber = 300;
particleTracking.LengthScalePe = 30;
particleTracking.VelocityScalePe = flowDFN.MeanVelocity;
particleTracking.MolecularDiffusion = particleTracking.LengthScalePe /
                                      particleTracking.PecletNumber *
                                      particleTracking.VelocityScalePe;
particleTracking.FactorTimeScaleCrossElement = 2;
particleTracking.TimeScaleCrossElement =
    pow(meshGen.MeanGridSize, 0.5) / flowDFN.MaxVelocity;
particleTracking.DeltaT = particleTracking.TimeScaleCrossElement /
                          particleTracking.FactorTimeScaleCrossElement;
particleTracking.FluexWeightedOrUniformInjection = true;
particleTracking.OutputAllPTInformationOrFPTCurve = true;
particleTracking.SpacingOfControlPlanes = 30;
particleTracking.IfOutputVarianceOfDisplacementsEachStep = true;
particleTracking.IfInjectAtCustomedPlane = true;
particleTracking.CustomedPlaneInjection = 23;
particleTracking.IfUseFluxWeightedOrEqualProbableMixingIntersection = true;
particleTracking.ParticleTracking(my_dfn, meshGen, flowDFN);
particleTracking.Visualization(my_dfn, meshGen, flowDFN,
                               "DFN_DISPERSION_VISUAL",
                               "DFN_DISPERSION_VISUAL", "DFN_FLOW_VISUAL");
```
An empty class of `cuDFNsys::PTDFN<double>` named `particleTracking` is generated firstly. The number of particles `particleTracking.NumParticles` is 20000. The number of time step `particleTracking.NumTimeSteps` is 200. The value of $Pe$ `particleTracking.PecletNumber` is 300. 

The length scale here `particleTracking.LengthScalePe` is the size of the domain in the $x$ direction (one can also define it as the expection of fracture sizes $R$). The velocity scale `particleTracking.VelocityScalePe` is the mean fracture velocity. The $D_m$ value `particleTracking.MolecularDiffusion` is then calculated based on these scales.

The time scale for a random walker to cross a finite element `particleTracking.TimeScaleCrossElement` is the square root of the mean size of finite elements over the maximum fracture velocity. The $\delta t$ value `particleTracking.DeltaT` is set to be `particleTracking.TimeScaleCrossElement` over a reduction factor `particleTracking.FactorTimeScaleCrossElement`.

The injection method `particleTracking.FluexWeightedOrUniformInjection = true` is flux-weighted. The spatial information of particles for each time step is wanted, and hence `particleTracking.OutputAllPTInformationOrFPTCurve = true.` Also, the variace of displacements of particles for each step is output by setting `particleTracking.IfOutputVarianceOfDisplacementsEachStep = true`.

The spacing of control planes is 30 by setting `particleTracking.SpacingOfControlPlanes = 30`.

Because molecular diffusion exists, the injection plane is better to be not at inlet, otherwise, particles will leave the domain from inlet. Therefore `particleTracking.IfInjectAtCustomedPlane = true` and `particleTracking.CustomedPlaneInjection = 23`, meaning that the injection plane is at $z = 23$ as the percolation direction is along the $z$ axis.

Finally, for the mixing rule at intersections, we set outgoing-flux-weighted mixing rule, so `particleTracking.IfUseFluxWeightedOrEqualProbableMixingIntersection = true`.

After the set-up of these member variables, we can run particle tracking by `particleTracking.ParticleTracking(my_dfn, meshGen, flowDFN)`. Note that if one wants to run more time steps in the same DFN, it is possible. The codes `QuickStartGuide_DFN_I_DFN.cu`, `QuickStartGuide_DFN_II_MESH.cu`, `QuickStartGuide_DFN_III_FLOW.cu` and `QuickStartGuide_DFN_IV_PT.cu` show us that each class can be implemented and stored in HDF5 files, and then these classes can be recovered by loading data in these HDF5 files. For example, to run more time steps, one just run the executable files `QuickStartGuide_DFN_I_DFN`, `QuickStartGuide_DFN_II_MESH` and `QuickStartGuide_DFN_III_FLOW` in order, then, run `QuickStartGuide_DFN_IV_PT` repeatly to obtain particle tracking results for longer time. 

The visualization scripts are output by `particleTracking.Visualization()`. One can run `DFN_DISPERSION_VISUAL.py` or `DFN_DISPERSION_VISUAL.m` to see the animation. But before that, the particle data should be addressed by an executable file, to transform the data to 3D. Just run `./Transform2DH5ParticleDataTo3D 0 DFN_MESH_VISUAL.h5` in ubuntu terminal. Then run `DFN_DISPERSION_VISUAL.py` or `DFN_DISPERSION_VISUAL.m`. 

One example is shown in the figure below.

<p align="center">
  <img width="500" src="https://github.com/qq1012510777/cuDFNsys/blob/main/Manual/A_DFN_PT.png">
</p>
<p align="center">
    <em>Visualization of a DFN particle tracking result (by matlab). The color and color bar denote the $z$ value of finite elements. </em>
</p>

#### 2.4.4. Particle tracking results (HDF5)
The results are stored in a directory named `ParticlePositionResult` under the current working directory. 

The file `DispersionInfo.h5` record the input parameters like molecular diffusion coefficient, $\delta t$, and so on.