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

In this code, a stochastic DFN with one group of fractures is generated. The number of fractures is `500`, the orientations of fracture follow the Fisher distribution, but the parameter $\kappa = 0$, meaning that the distribution is actually uniform. The sizes of fractures follow a power-law distribution, with $\alpha = 1.5$, the minimum and maximum sizes of fractures are 1 and 15, respectively. Note that the size of fractures is described by the radius $R$ of circumscribed circles of the square fractures.

The condutivity of fractures is related to $R$. The aperture $b$ of fractures is 
$$b = \gamma R ^\beta,$$
where $\gamma $ and $\beta$ are a constant and an exponent, respectively. The conductivity $k_f$ of fractures is
$$k_f = \frac{b^3}{12}$$.
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

The variable `ParaSizeDistri` (type: `cuDFNsys::Vector4<double>`) is a vector of four elements. It is equivalent to `double4` or `float4` in cuda, depending on the template you set. Here the template is `<double>`, so `ParaSizeDistri` is actually a `double4` with elements that can be accessed by `ParaSizeDistri.x`, `ParaSizeDistri.y`, `ParaSizeDistri.z` and `ParaSizeDistri.w`. To create a `cuDFNsys::Vector4<double>`, we can use a _cuDFNsys_ function: `cuDFNsys::MakeVector4(double, double, double, double)`.

When `ModeSizeDistri = 0`, `ParaSizeDistri.x` is $\alpha$, `ParaSizeDistri.y` is the minimum $R$, `ParaSizeDistri.z` is the maximum $R$, `ParaSizeDistri.w` means nothing.

