# cuDFNsys

_cuDFNsys_ is an open-source CUDA library (under the GPL license) for DFN generations. It is also able to simulate flow (steady-state) and particle transport in DFNs, based on the mixed hybrid finite element method and particle tracking algorithm. _cuDFNsys_ contains around 10 000 lines of codes (counted on June 29th, 2023).

_cuDFNsys_ does not support the GUI and can only run on Ubuntu. _cuDFNsys_ is not friendly to users unfamiliar with Linux, C++, Make and CMake.

<p align="center">
  <img width="500" src="https://github.com/qq1012510777/cuDFNsys/blob/main/moive_particlesII.gif">
</p>
<p align="center">
    <em>Figure 1. Dispersion in a DFN with 300 000 particles. The movement of particles is purely driven by the advection. Left: the DFN and hydraulic head field. Right: the particle spreading. </em>
</p>

<p align="center">
  <img width="500" src="https://github.com/qq1012510777/cuDFNsys/blob/main/moive_particles_Diffusion_dominated.gif">
</p>
<p align="center">
    <em>Figure 2. Diffusion-dominated dispersion in a DFN with 300 000 particles. The Peclet number is 0.01. The length scale in the Peclet number is the mean value of fracture sizes. Left: the DFN and hydraulic head field. Right: the particle spreading. The particles initialize at a plane with z = 0. </em>
</p>

# Authors

* Tingchang YIN, Westlake University & Zhejiang University, China, yintingchang@foxmail.com

* Sergio GALINDO-TORRES, Westlake University, China, s.torres@westlake.edu.cn

# Terms

_cuDFNsys_ - simulating flow and transport in 3D fracture networks
Copyright (C) 2022, Tingchang YIN, Sergio GALINDO-TORRES 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Prerequisites and installation
_cuDFNsys_ should be installed and run on Ubuntu.

_cuDFNsys_ relies on several open packages: OpenMP, CUDA, Eigen, Gmsh, UMFPACK and HDF5.

The Gmsh C++ API that _cuDFNsys_ relies on should support the OCC mode, meaning that the occt library is required.

The following is the Ubuntu command lines to install the relying packages:
```
# first update the package repository information
cd ~
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install cmake 
# note that the version of cmake should be 3.10 or higher
# if the version is not satisfied, try to install cmake from source

# cuda
sudo apt-get install nvidia-cuda-dev
sudo apt-get install nvidia-cuda-toolkit 

# openMP
sudo apt-get install libomp-dev

# install Eigen
sudo apt-get install libeigen3-dev
# create link (optional)
sudo ln -s /usr/include/eigen3/Eigen/  /usr/include/Eigen

# install blas, lapack and umfpack (suitesparse)
sudo apt-get install libblas-dev liblapack-dev
sudo apt-get install libsuitesparse-dev

# install HDF5
sudo apt-get install libhdf5-dev

# install occt for gmsh
sudo apt-get install libocct-foundation-dev libocct-data-exchange-dev

# install fltk for gmsh
sudo apt-get install libfltk1.3-dev

# now install gmsh from source
git clone http://gitlab.onelab.info/gmsh/gmsh.git
cd gmsh
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=path-to-install-gmsh  ..
make
sudo make install  

# change the path-to-install-gmsh to where you want to install the gmsh

```
After the installation of these relying packages (by `sudo apt-get install`, they are install at default locations), the _cuDFNsys_ package can be installed by the following steps:
```
cd ~/cuDFNsys/lib
# here I put the cuDFNsys source code under HOME
```
Open `CMakeLists.txt`, find Line 28: `set (CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")`, change `"/usr/local/cuda/bin/nvcc"` to the path of `nvcc` in your computer.

Then, open `SetPaths.cmake`, you will see
```
# Eigen include directory
SET(EIGEN_INCLUDE_SEARCH_PATH   $ENV{HOME}/somewhere)

# Gmsh include directory
SET(GMSH_INCLUDE_SEARCH_PATH    $ENV{HOME}/somewhere)

# Gmsh lib directory
SET(GMSH_LIBRARY_SEARCH_PATH    $ENV{HOME}/somewhere)

# umfpack include directory
SET(UMFPACK_INCLUDE_SEARCH_PATH $ENV{HOME}/somewhere)

# umfpack lib directory
SET(UMFPACK_LIBRARY_SEARCH_PATH $ENV{HOME}/somewhere)
```
Again, change these paths like `$ENV{HOME}/somewhere` to the path of the packages in your computer, e.g., `/usr/include/eigen3`.

Now, you can compile the _cuDFNsys_:
```
cd ~/cuDFNsys/lib
bash compileCode.sh
```
If this script is finished without errors, then the file `libcuDFNsys.a` appears in `~/cuDFNsys/lib`.

More details about the installation of the relying packages:

* [Installation of CUDA](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html)

* [Installation of Gmsh and occt](https://gitlab.onelab.info/gmsh/gmsh/-/wikis/Gmsh-compilation)

* [Installation of Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download)

* [Installation of UMFPACK](https://github.com/DrTimothyAldenDavis/SuiteSparse)

# Compile a quickstart example
By `cd ~/cuDFNsys/QuickStartGuide`, a `Makefile` can be seen there. Open it, you can see
```
# NVCC path
NVCC=/usr/local/cuda/bin/nvcc
# include paths for headers
Hdf5IncludePath=/usr/lib/x86_64-linux-gnu/hdf5/serial/include
GmshIncludePath=$(HOME)/pkg/gmsh-4.8.4-source/MY_GMSH/include
EigenIncludePath=$(HOME)/pkg/eigen
UmfpackIncludePath=$(HOME)/pkg/SuiteSparse-master/include
# library paths
GmshLibraryPath=$(HOME)/pkg/gmsh-4.8.4-source/MY_GMSH/lib
UmfpackLibraryPath=$(HOME)/pkg/SuiteSparse-master/lib
Hdf5LibraryPath=/usr/lib/x86_64-linux-gnu/hdf5/serial
```
Change the paths of the headers and libraries, and compile the QuickStartGuide.cu just by
```
make
```

# Visualization

After simulation, _cuDFNsys_ outputs .h5 and .m (and/or .py) files, and you can visualize them by the generated .m or .py file. 

The visulization with MATLAB is simple, you just need the license of MATLAB. 

With Python, you need to install the _mayavi_ Python package. _mayavi_ is a very good visualization engine. The Matplotlib package is not good at 3D visualizations, see [Why my 3D plot doesn’t look right at certain viewing angles](https://matplotlib.org/2.2.2/mpl_toolkits/mplot3d/faq.html).

# Manual
Right now, I provide a quickstart guide to explain how can one do simulations with _cuDFNsys_ functions. See QuickStartGuide/src/QuickStartGuide.cu.

To compile and run QuickStartGuide example, run `make` in the directory-'QuickStartGuide'.

Then, run the command in your shell: `./QuickStartGuide 10000`

The number `10000` above means that 10000 particles are injected.

Note that the particles' information (e.g., positions) are 2D. You have to transform them into 3D data with the command `./Transform2DH5ParticleDataTo3D 0 DFN_mesh.h5`.

# Directories

_QuickStartGuide_: a quickstart guide CUDA example to show how to do simulations with _cuDFNsys_ functions.

_Modules_: cmake script to find packages that _cuDFNsys_ relies on in a Ubuntu system

_include_: header files containing declarations of functions/classes

_src_: source files containing definitions of functions/classes

Other directories are all about my personal examples (user's interfaces), with different purposes, meanwhile showing that how to call _cuDFNsys_ functions, with `make` and other tools.