# cuDFNsys

The _cuDFNsys_ is an open-source C++/CUDA library (under GPL license) for DFN generations, based on CUDA. It is also able to simulate static flow and particle transport in DFN, based on mixed hybrid finite element method and particle tracking algorithm.

_cuDFNsys_ does not support GUI and can only run on Ubuntu. _cuDFNsys_ is not friendly to users who do not familiar with Linux, C++ and CMake.

_cuDFNsys_ right now can only generate one group of fractures, because I am studying percolation in DFNs and some quantities in percolation theory is easier to calculate with only one fracture group. Geometrical attributes of fractures follow certain distributions. In reality, multiple fracture famlies might be observed.

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

# Prerequisites
The _cuDFNsys_ should be installed and run on Ubuntu.

The _cuDFNsys_ relies on several open packages: CUDA, Eigen, Gmsh, UMFPACK and HDF5.

Installations of these libraries could be disturbing, even you are familiar with Linux, CMake, C++ and so on. I cannot write a bash file that help you install all these libraries at one time. You should go to their homepages, and install them one by one. You may feel painful when install some libraries, for instance, the cuda package is difficult to install, and problem might happen occasionally. The Gmsh C++ API that _cuDFNsys_ relies on should support OCC mode, meaning that the occt library is required. Also, the installation of UMFPACK could be difficult. _Anyway, I am very willing to help you install them, but you should be familiar with Linux, CMake, C++ and so on_.

* [Installation of CUDA](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html)

* [Installation of Gmsh and occt](https://gitlab.onelab.info/gmsh/gmsh/-/wikis/Gmsh-compilation)

* [Installation of Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download)

* [Installation of UMFPACK](https://github.com/DrTimothyAldenDavis/SuiteSparse)

# Environment setup
After installation of above-mentioned packages, the _cuDFNsys_ should be linked to them. In directory 'Modules', .cmake files can be found, where the path of corresponding package should be added/set.

For example, in 'Modules/FindGMSH.cmake', you'll see

    SET(GMSH_INCLUDE_SEARCH_PATH
        $ENV{HOME}/pkg/gmsh-4.8.4-source/MY_GMSH/include
    )
    
    SET(GMSH_LIBRARY_SEARCH_PATH
        $ENV{HOME}/pkg/gmsh-4.8.4-source/MY_GMSH/lib
    )

The two user-set variables are paths of Gmsh header file and libraries. You have to change them if you install Gmsh in other directory.

In 'QuickStartGuide/CMakeLists.txt', you'll see

    set (CUDFNSYS_ROOT $ENV{HOME}/cuDFNsys)

The above variable is the path of _cuDFNsys_. You, again, might change it.

# Visualization

After simulation, _cuDFNsys_ outputs .h5 files, and you can visualize them by the generated .m file or .py file. 

Visulization with MATLAB is simple, you just need the license of MATLAB. 

With Python, you need to install the mayavi Python package. Sorry again for the inconvenience, the installation of mayavi could be difficult, but the mayavi is a very good visualization engine. The Matplotlib package is not good at 3D visualization, see [Why my 3D plot doesn’t look right at certain viewing angles](https://matplotlib.org/2.2.2/mpl_toolkits/mplot3d/faq.html).

# Manual
Right now, I provide a quickstart guide to explain how can one do simulation with _cuDFNsys_ functions. See QuickStartGuide/src/main.cu.

To compile and run QuickStartGuide example, run compileCode.sh in the directory-'QuickStartGuide'.

For example, run the command with your shell: **bash compileCode.sh 1000 500 1e5 0**

Explanation of the four arguments:

* Number of particles: 1000

* Number of time steps: 500

* Delta t: 1e5

* Diffusion (local): 0

# Directories

_QuickStartGuide_: a quickstart guide CUDA example to show how to do simulation with _cuDFNsys_ functions.

_CodingGuideline_: my coding guidline, e.g. how to name a member variable in a class

_Modules_: cmake script to find packages that the _cuDFNsys_ relies on in a Ubuntu system

_include_: header files containing declarations of functions/classes

_src_: source files containing definitions of functions/classes

_PercolationTest_, _TestFunctionsOnlyForDevlopment_, _TestResolutionEffect_ and other directories are all about examples (user's interfaces), with different purposes, meanwhile showing that how to call _cuDFNsys_ functions, with cmake and other tools.