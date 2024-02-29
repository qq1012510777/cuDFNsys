# cuDFNsys

_cuDFNsys_ is an open-source CUDA library (under the GPL license) for DFN generations. It is also able to simulate flow (steady-state) and particle transport in DFNs, based on the mixed hybrid finite element method and particle tracking algorithm. _cuDFNsys_ contains around 10 000 lines of codes (counted on June 29th, 2023).

_cuDFNsys_ does not support the GUI and can only run on Ubuntu.

_cuDFNsys_ is an easy-implemented, object-oriented CUDA C++ library. The simulation is performed just by establishing CUDA C++ classes, defining member variables and calling member functions. The results are stored in a HDF5 format. At runtime, the fracture, mesh or flow data can be accessed by visting member variables of classes.

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

* [Tingchang YIN](https://qq1012510777.github.io/), Westlake University & Zhejiang University, China, yintingchang@foxmail.com

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

# Relevant publications
- [T Yin, T Man, Ling Li, SA Galindo-Torres. Finite-Size Scaling for the Permeability of Discrete Fracture Networks (2023). Geophysical Research Letters, 50, e2022GL100837](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022GL100837)
- [T Yin, T Man, SA Galindo-Torres. Universal scaling solution for the connectivity of discrete fracture networks (2022). Physica A: Statistical Mechanics and its Applications 599, 127495](https://www.sciencedirect.com/science/article/abs/pii/S0378437122003557)

# Prerequisites
_cuDFNsys_ should be installed and run on Ubuntu.

_cuDFNsys_ relies on several open packages: OpenMP, CUDA, Eigen, Gmsh, UMFPACK and HDF5. These dependencies can be installed by using `sudo apt-get install packageName-dev`.

The following is the Ubuntu command lines to install the relying packages:
```
# run the following command one by one to install dependencies
# I have tested them on Ubuntu 23.04
# first update the package repository information
cd ~
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install gfortran
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

# now install gmsh
sudo apt-get install libgmsh-dev
```

# Building cuDFNsys

After the installation of dependencies (by `sudo apt-get install`, these relying packages are installed at default locations), and the _cuDFNsys_ library can be installed by the following steps:
```
git clone https://github.com/qq1012510777/cuDFNsys.git
cd ~/cuDFNsys/lib
# here I put the cuDFNsys source code under HOME
```
Open `SetPaths.cmake`, one will see
```
# NVCC directory
SET(NVCC_PATH                   /usr/lib/nvidia-cuda-toolkit/bin/nvcc)

# Eigen include directory
SET(EIGEN_INCLUDE_SEARCH_PATH   /usr/include)

# Gmsh include directory
SET(GMSH_INCLUDE_SEARCH_PATH    /usr/include)

# Gmsh lib directory
SET(GMSH_LIBRARY_SEARCH_PATH    /usr/lib/x86_64-linux-gnu)

# umfpack include directory
SET(UMFPACK_INCLUDE_SEARCH_PATH /usr/include/suitesparse)

# umfpack lib directory
SET(UMFPACK_LIBRARY_SEARCH_PATH /usr/lib/x86_64-linux-gnu)
```
The above script of CMake shows that paths of the dependencies in a computer with Ubuntu 23.04. Generally these paths do not need to be changed.

But, if one install the dependencies in different paths, change these paths, e.g., `/usr/include` to the path of the packages in your computer, e.g., `/path/in/your/computer`.

Now, one can compile _cuDFNsys_ by:
```
cd ~/cuDFNsys/lib
mkdir build 
cd build 
cmake -DCMAKE_CUDA_ARCHITECTURES=60  ..
make
cd ..
```
If these commands are finished without errors, then the file `libcuDFNsys.a` should appear in `~/cuDFNsys/lib`.

More details about the installation of the relying packages:

* [Installation of CUDA](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html)

* [Installation of Gmsh and occt](https://gitlab.onelab.info/gmsh/gmsh/-/wikis/Gmsh-compilation)

* [Installation of Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download)

* [Installation of UMFPACK](https://github.com/DrTimothyAldenDavis/SuiteSparse)

# Compile a quickstart example
By `cd ~/cuDFNsys/QuickStartGuide`, a `Makefile` can be seen there. Open it, and one can see
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
Change the paths of the headers and libraries (if necessary), and compile the QuickStartGuide.cu just by
```
make
```
and seversal executable files can be generated.

After the compilation of QuickStartGuide, run one executable file by
```
./QuickStartGuide
```
If errors happen, e.g., `error while loading shared libraries: libgmsh.so: cannot open shared object file: No such file or directory`, just add two enviromental variables to `~/.bashrc` as follows:
```
vim ~/.bashrc
```
add the following content to `~/.bashrc`
```
export LD_LIBRARY_PATH=path-to-gmsh-library:$LD_LIBRARY_PATH
export LIBRARY_PATH=path-to-gmsh-library:$LIBRARY_PATH
# change path-to-gmsh-library to the dynamic gmsh library in the computer
```
then update `~/.bashrc` by
```
source ~/.bashrc
```
This quickstart example can be run by `./QuickStartGuide`, where a DFN is generated in which the flow is solved by a mixed hybrid finite element method and the particle tracking is implement for a not-very-long time. 

_cuDFNsys_ runs the flow and transport in DFNs by four cuda c++ classes, and the DFN model and mesh and flow data are reproducable, for example, to run more steps for particle tracking in one DFN, please see codes (i.e., `ReRunPT.cpp``) and run the corresponding executable file. More details are in [manual](Manual/Manual.md).

# Visualization

After simulation, _cuDFNsys_ outputs .h5 and .m (and/or .py) files, and you can visualize them by the generated .m or .py file. 

The visulization with MATLAB is simple, you just need a MATLAB installed, and run the `.m` file. 

With Python, you need to install the _mayavi_ Python package. _mayavi_ is a very good visualization engine. The Matplotlib package is not good at 3D visualizations, see [Why my 3D plot doesn’t look right at certain viewing angles](https://matplotlib.org/2.2.2/mpl_toolkits/mplot3d/faq.html).

Installation of _mayavi_ can be done by the following commands (I've test them in Ubuntu 23.04):
```
sudo apt install python3-vtk9 python3-pip python3-pyqt5 python3-traits python3-traitsui python3-numpy python3-matplotlib python3-setuptools python3-pyqt5.qtopengl
# Please note that package names and dependencies may change over time, so make sure to check for any updates or changes in package names specific to your Ubuntu version.
sudo apt-get install mayavi2
```
To run python visualization script, just run `python3 NameOfScript.py`.

# Manual
Manual for _cuDFNsys_ is [here](Manual/Manual.md).

Now, _cuDFNsys_ has a simple GUI, which rquires the _cuDFNsys_ library to be compiled. Also, the GUI can visualize the DFN, mesh, flow, and PT both in Python Matplotlib and Mayavi (required it to be installed). See [Manual](Manual/Manual.md) for more details about how to run a cuDFNsys GUI.

# Directories

_QuickStartGuide_: a quickstart guide CUDA example to show how to do simulations with _cuDFNsys_ functions.

_Modules_: cmake script to find packages that _cuDFNsys_ relies on in a Ubuntu system

_include_: header files containing declarations of functions/classes

_src_: source files containing definitions of functions/classes

_SomeApplications_: these are all about my personal examples (user's interfaces), with different purposes, meanwhile showing that how to call _cuDFNsys_ functions, with `make` and other tools.