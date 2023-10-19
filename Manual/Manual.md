# Manual for _cuDFNsys_

* Name: Tingchang YIN

* Institution: Westlake University, China

* Date: Oct. 19, 2023

* Update date: Oct. 19, 2023 

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
The above code just generates a vector (type: `cuDFNsys::Vector4<double>`) named `vec4` of four elements (double precision), which can be accessed by `vec4.x`, `vec4.y`, `vec4.z` and `vec4.w`. More details about `cuDFNsys::Vector4<double>` are discussed later.

To compile this code with `libcuDFNsys.a`, we can copy the following `Make` script to `~/cuDFNsys/CompilationTest/Makefile`:
```
ExeName=CompilationTest
NVCC=/usr/local/cuda/bin/nvcc
# include paths for headers
cuDFNsysIncludePath=$(HOME)/cuDFNsys/include
Hdf5IncludePath=/usr/lib/x86_64-linux-gnu/hdf5/serial/include
GmshIncludePath=$(HOME)/pkg/gmsh-4.8.4-source/MY_GMSH/include
EigenIncludePath=$(HOME)/pkg/eigen
UmfpackIncludePath=$(HOME)/pkg/SuiteSparse-master/include
# library paths
cuDFNsysLibraryPath=$(HOME)/cuDFNsys/lib
GmshLibraryPath=$(HOME)/pkg/gmsh-4.8.4-source/MY_GMSH/lib
UmfpackLibraryPath=$(HOME)/pkg/SuiteSparse-master/lib
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
NVCC=/usr/local/cuda/bin/nvcc
# include paths for headers
cuDFNsysIncludePath=$(HOME)/cuDFNsys/include
Hdf5IncludePath=/usr/lib/x86_64-linux-gnu/hdf5/serial/include
GmshIncludePath=$(HOME)/pkg/gmsh-4.8.4-source/MY_GMSH/include
EigenIncludePath=$(HOME)/pkg/eigen
UmfpackIncludePath=$(HOME)/pkg/SuiteSparse-master/include
# library paths
cuDFNsysLibraryPath=$(HOME)/cuDFNsys/lib
GmshLibraryPath=$(HOME)/pkg/gmsh-4.8.4-source/MY_GMSH/lib
UmfpackLibraryPath=$(HOME)/pkg/SuiteSparse-master/lib
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

If errors happen when you run `./main`, e.g., `error while loading shared libraries: libgmsh.so: cannot open shared object file: No such file or directory`, we just add two enviromental variables to `~/.bashrc` as follows:
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

## 2. Documented _cuDFNsys_ functions/classes 
Here I document some _cuDFNsys_ functions/classes (which are directly related to the generation of DFNs, flow and transport simulations). I'll do it with an introduction of a quickstart example, which is in `~/cuDFNsys/QuickStartGuide/src/QuickStartGuide.cu`.

I will explain this code line by line, accompanying by the description of _cuDFNsys_ functions/classes.

## 2.1 DFN generation
In this code, a stochastic DFN with one group of fractures is generated. The number of fractures is `500`, the orientations of fracture follow the Fisher distribution, but the parameter $\kappa = 0$, meaning that the distribution is actually uniform. The sizes of fractures follow a power-law distribution, with $\alpha = 1.5$, the minimum and maximum sizes of fractures are 1 and 15, respectively. Note that the size of fractures is described by the radius $R$ of circumscribed circles of the square fractures.

The condutivity of fractures is related to $R$. The aperture $b$ of fractures is 
$$b = \gamma R ^\beta,$$
where $\gamma$ and $\beta$ are a constant and an exponent, respectively. The conductivity $k_f$ of fractures is
$$k_f = \frac{b^3}{12}.$$
Therefore, the conductivity of fractures depends on the size of fractures.

We set $\gamma = 5.5e-4, \beta = 0.25$. Therefore, the following variables describe the fracture parameters

```
    int NumFractures = 500;
    double kappa = 0;
    double beta = 0.25, gamma = 5.5e-4;
    int ModeSizeDistri = 0;
    cuDFNsys::Vector4<double> ParaSizeDistri =
        cuDFNsys::MakeVector4(1.5,
                              1.,
                              15.,
                              0.);
```
The variable `ModeSizeDistri` (an integer) is an indicator for the size distributions: 0 for power-law, 1 for lognormal, 2 for uniform, 3 for mono-sized. 

The variable `ParaSizeDistri` (type: `cuDFNsys::Vector4<double>`) is a vector of four elements. It is equivalent to `double4` or `float4` in cuda, depending on the template you set. That is to say, the template can be `double` or `float`. Here the template is `<double>`, so `ParaSizeDistri` is actually a `double4` with elements that can be accessed by `ParaSizeDistri.x`, `ParaSizeDistri.y`, `ParaSizeDistri.z` and `ParaSizeDistri.w`. To create a `cuDFNsys::Vector4<double>`, we can use a _cuDFNsys_ function: `cuDFNsys::MakeVector4(double, double, double, double)`.

When `ModeSizeDistri = 0`, the size distribution is a power-law. `ParaSizeDistri.x` is $\alpha$, `ParaSizeDistri.y` is the minimum $R$, `ParaSizeDistri.z` is the maximum $R$, `ParaSizeDistri.w` means nothing.

When `ModeSizeDistri = 1`, the size distribution is lognormal `ParaSizeDistri.x` is the mean of logarithmic values, `ParaSizeDistri.y` is the standard deviation of logarithmic values, `ParaSizeDistri.z` is the minimum $R$, `ParaSizeDistri.w` is the maximum $R$.

When `ModeSizeDistri = 2`, the size distribution is uniform. `ParaSizeDistri.x` is the minimum $R$, `ParaSizeDistri.y` is the maximum $R$, `ParaSizeDistri.z` and `ParaSizeDistri.w` mean nothing.

When `ModeSizeDistri = 3`, the size distribution is mono-sized. `ParaSizeDistri.x` is $R$, `ParaSizeDistri.y`, `ParaSizeDistri.z` and `ParaSizeDistri.w` mean nothing.

Then, I want the domain be a cuboid column. The domain size in x is 30, size in y is 30, size in z is 60. The parameters for domain are as follows:
```
    double DomainSize_X = 30;
    double3 DomainDimensionRatio = make_double3(1, 1, 2);
    int perco_dir = 2;
```
`DomainDimensionRatio` is a `double3` variable (that is a built-in data-type in cuda, also the function `make_double3` is a built-in function in cuda), and it clarifies the dimension ratio with respect to the size in x. Here the ratio is `(1, 30/30 = 1, 60 / 30 = 2)`. `perco_dir` defines the percolation direction, i.e., the mean flow direction. We later check the percolation state in the z directon. So, `perco_dir = 0` means the x direction, `perco_dir = 1` means the y direction, and `perco_dir = 2` means the z direction.

Next, I create an empty vector of structs (the fracture) by
```
thrust::host_vector<cuDFNsys::Fracture<double>> Frac_verts_host(NumFractures);
```
`thrust::host_vector<type>` is a built-in container in cuda, which creates a vector of something on the host side (manipulations on CPU). `cuDFNsys::Fracture<double>` is a _cuDFNsys_ struct, with templates (`double` or `float`). A `cuDFNsys::Fracture<double>` has several member variables:
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
`cuDFNsys::Vector3<T>` is `double3` or `float3`, depending on the template you set. `Conductivity` is $k_f$. `Verts3D` is the four vertices of the square fractures. `center` is the center of the fracture. `Verts3DTruncated` is the vertices of square fractures after the fracture is truncated by the domain. `NumVertsTruncated` is the number of vertices after truncation, so its maximum number is eight.  `Radius` is the radius of fractures. `ConnectModelSurf` denotes if the fracture intersects the six surfaces of the domain, `true` means intersected. `NormalVec` is a `double3` or `float3`, denoting the normal vector (always pointing up) of a fracture. 

The fracture vector is named `Frac_verts_host` with size of `NumFractures`. To generate fractures on GPU side, I create a device vector and a pointer to this device vector.
```
    thrust::device_vector<cuDFNsys::Fracture<double>> Frac_verts_device(NumFractures);
    cuDFNsys::Fracture<double> *Frac_verts_device_ptr;
    Frac_verts_device_ptr = thrust::raw_pointer_cast(Frac_verts_device.data());
```
`thrust::raw_pointer_cast` is a built-in function in the `thrust` function in cuda. It links the address of a thrust vector to a pointer. For more details, one can google it.

Now, we have all parameter defined in different variables. Let's generate a stochastic DFN.
```
    time_t t;
    time(&t);
    cuDFNsys::Fractures<double><<<NumFractures / 256 + 1, 256>>>(Frac_verts_device_ptr,
                                                                 (unsigned long)t,
                                                                 NumFractures,
                                                                 DomainSize_X,
                                                                 ModeSizeDistri,
                                                                 ParaSizeDistri,
                                                                 kappa,
                                                                 beta,
                                                                 gamma,
                                                                 DomainDimensionRatio);

    cudaDeviceSynchronize();
    Frac_verts_host = Frac_verts_device;
```
`time_t t;` and `time(&t);` are used to address random seeds. `t` is now the random seed. `cuDFNsys::Fractures<double>` is a global function to generate fractures on GPU, with template of `double`. `<<<NumFractures / 256 + 1, 256>>>` is a configuration of the execution parameters for a CUDA kernel launch:
```
template <typename T>
__global__ void Fractures(cuDFNsys::Fracture<T> *verts,
                          unsigned long seed,
                          int count,
                          T model_L,
                          uint ModeSizeDistri,                 
                          cuDFNsys::Vector4<T> ParaSizeDistri,
                          T kappa,
                          T beta,
                          T gamma,
                          double3 DomainDimensionRatio = make_double3(1, 1, 1),
                          cuDFNsys::Vector3<T> MeanOrientation = cuDFNsys::MakeVector3((T)0., (T)0., (T)1.));
```
* `*verts` is the pointer to the device fracture vector.
* `seed` is a random seed, with type casting to `unsigned long`
* `count` is a integer, which is the number of fractures
* `model_L` is the domain size in x direction
* `ModeSizeDistri` is a integer, which is an indicator for fractur size distributions
* `ParaSizeDistri` is a `cuDFNsys::Vector4<double>` here, meaning the parameters for the size distribution
* `kappa` is a double, meaning the $\kappa$ value for the Fisher distribution
* `beta` is a double, meaning the $\beta$ value for conductivity
* `gamma` is a double, meaning the $\gamma$ value for conductivity
* `DomainDimensionRatio` is a `cuDFNsys::Vector3<double>` here, denoting the dimension ratio of the domain size, the default value is `(1, 1, 1)`
* `MeanOrientation` is the mean orientation of the Fisher distribution, the defauly value is `(0, 0, 1)`

`cudaDeviceSynchronize()` is a CUDA runtime API function used in CUDA C/C++ programming to ensure that all previously issued CUDA kernel launches have completed and that the GPU device has synchronized with the host CPU. In other words, it's used to make sure that all GPU work is finished before proceeding with further CPU operations.

After the kernel function, we copy data from device to host by `Frac_verts_host = Frac_verts_device;`. Now the generate of fracture networks is finished.

## 2.2 Identification of fracture intersections
To identify the intersetion of fractures, we firstly create an empty map:
```
    std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<double>, cuDFNsys::Vector3<double>>> Intersection_map;
```
The key of `Intersection_map` is a pair of `size_t`, which is the pair of IDs of fractures that intersect with each other. The value of the map is also a pair, which is a pair of `cuDFNsys::Vector3<double>`, denoting the coordinates of two ends of the intersections. The next _cuDFNsys_ class uses this map as an input parameter
```
    cuDFNsys::IdentifyIntersection<double> identifyInters{Frac_verts_host.size(),
                                                          Frac_verts_device_ptr,
                                                          true,
                                                          Intersection_map};
```
The constructor of this class needs four parameters:
```
template <typename T>
class IdentifyIntersection
{
public:
    // constructor GPU
    IdentifyIntersection(const size_t &NumOfFractures,
                         cuDFNsys::Fracture<T> *Frac_verts_device_ptr,
                         const bool &If_Trucncated,
                         std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> &Intersection_map);
};
```
* `NumOfFractures` is the number of fractures
* `Frac_verts_device_ptr` is a pointer to the device vector of fractures
* `If_Trucncated` means that if the trucated fractures or original fractures are considered in identification of intersections
* `Intersection_map` is the map to store intersection pairs and intersection coordinates.

## 2.3 Identification of fracture clusters
A cluster means that a number of fractures connect each other one by one. If no other fractures connect to one fracture, then this fracture is also a cluster.

We can identify the cluster in DFNs, and check if there is at least one cluster spanning the domain in the pre-defined percolation direction, i.e., percolation clusters.

Two empty vectors are created:
```
    std::vector<std::vector<size_t>> ListClusters;
    std::vector<size_t> Percolation_cluster;
```
Both are `std::vector`. The first one is `std::vector` of `std::vector`. It will store the clusters in DFNs. For instance, `ListClusters[0]`, `ListClusters[1]` and so on are clusters in a DFN. `ListClusters[0][0]`, `ListClusters[0][1]` and so on are IDs of fractures that belong to the same cluster.

`Percolation_cluster` will store the ID of the CLUSTERs that are percolative. For example, if `Percolation_cluster[0] = 0`, it means `ListClusters[0]` is a percolative cluster.

The identification is implement by the following class and functions
```
    cuDFNsys::Graph<double> G{(size_t)NumFractures, Intersection_map};
    G.UseDFS(ListClusters);
```
The constructor of class `cuDFNsys::Graph<double>` needs two parameters:
```
template <typename T>
class Graph
{
    // constructor
    Graph(const size_t &NumOfFractures,
          std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> Intersection_map);
};
```
* `NumOfFractures` is the number of fractures
* `Intersection_map` is the map to store intersection pairs and intersection coordinates
Then, the member function `UseDFS` needs one parameter:
```
template <typename T>
class Graph
{
    void UseDFS(std::vector<std::vector<size_t>> &ListClusters);
};
```
* `ListClusters` is the empty `std::vector` of `std::vector`. It will store the clusters in DFNs.

The percolation cluster can be identified by
```
    cuDFNsys::IdentifyPercolationCluster<double> IdentiClu{ListClusters,
                                                           Frac_verts_host,
                                                           perco_dir,
                                                           Percolation_cluster};
```
`IdentiClu` is a class of `cuDFNsys::IdentifyPercolationCluster<double>`:
```
template <typename T>
class IdentifyPercolationCluster
{
public:
    // constructor
    IdentifyPercolationCluster(const std::vector<std::vector<size_t>> &ListClusters,
                               const thrust::host_vector<cuDFNsys::Fracture<T>> &Frac_verts_host,
                               const int &Perco_dir,
                               std::vector<size_t> &Percolation_cluster);
};
```
* `ListClusters` stores the clusters in DFNs.
* `Frac_verts_host` is the host vector of fractures
* `Perco_dir` is the pre-defined percolation direction
* `Percolation_cluster` stores the IDs of percolative cluster

## 2.4 Visualization of DFNs
The visualization of DFNs can be implemented by:
```
    cuDFNsys::MatlabPlotDFN<double> PlotDFN_{"DFN.h5",
                                             "DFN.m",
                                             Frac_verts_host,
                                             Intersection_map,
                                             ListClusters,
                                             Percolation_cluster,
                                             false,
                                             true,
                                             true,
                                             true,
                                             DomainSize_X,
                                             perco_dir,
                                             true,
                                             "DFN",
                                             DomainDimensionRatio};
```
`PlotDFN_` is a class of `cuDFNsys::MatlabPlotDFN<double>`:
```
template <typename T>
class MatlabPlotDFN
{
public:
    // constructor
    MatlabPlotDFN(string h5_file_name,                                                                                    
                  string m_file_name,                                                                                 
                  thrust::host_vector<cuDFNsys::Fracture<T>> Frac_verts_host,                                         
                  std::map<pair<size_t, size_t>, pair<cuDFNsys::Vector3<T>, cuDFNsys::Vector3<T>>> Intersection_map, 
                  std::vector<std::vector<size_t>> ListClusters,                                                      
                  std::vector<size_t> Percolation_cluster,                                                            
                  bool If_show_truncated_frac,                                                                       
                  bool If_show_intersection,
                  bool If_show_cluster,
                  bool If_show_orientation,
                  T DomainSize_X,
                  int perco_dir,
                  bool if_python_visualization = false,
                  string PythonName_Without_suffix = "DFN_py",
                  double3 DomainDimensionRatio_d = make_double3(1, 1, 1));
};
```
* `h5_file_name` denotes the name of the `.h5` file, the suffix '.h5' is needed.
* `m_file_name` denotes the name of the `.m` file. It is a matlab script that will be generated in the working directory. The suffix '.m' is needed
* `Frac_verts_host` is the host vector of fractures
* `Intersection_map` is the map to store intersection pairs and intersection coordinates
* `ListClusters` stores the clusters in DFNs.
* `Percolation_cluster` stores the IDs of percolative cluster
* `If_show_truncated_frac` is a bool variable. `true` means displaying truncated fractures.
* `If_show_intersection` is a bool variable. `true` means displaying intersections.
* `If_show_cluster` is a bool variable. `true` means displaying clusters in different colors.
* `If_show_orientation` is a bool variable. `true` means displaying the distribution of fracture orientations.
* `DomainSize_X` is the domain size in x direction
* `perco_dir` is the pre-defined percolation direction. It should be `0`, `1` or `2`.
* `if_python_visualization` is a bool variable. `true` means outputting python script to visualize the DFN.
* `PythonName_Without_suffix` denotes the name of the python script, without the suffix `.py`
* `DomainDimensionRatio` is a `cuDFNsys::Vector3<double>` here, denoting the dimension ratio of the domain size, the default value is `(1, 1, 1)`

A matlab script named "DFN.m" will be generated in the current working directory. We can run it in Matlab. We can get the visualization of the DFN:
<p align="center">
  <img width="500" src="https://github.com/qq1012510777/cuDFNsys/blob/main/Manual/VISUALIZE_DFNs.png">
</p>
<p align="center">
    <em>Visualization of a DFN. From left to right: the DFN, the DFN with clusters in different colors, and the distribution of orientations. </em>
</p>
