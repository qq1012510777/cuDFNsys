# cuDFNsys

The _cuDFNsys_ is an open-source C++/CUDA library (under GPL license) for DFN generations, based on CUDA. It is also able to simulate static flow and particle transport in DFN, based on mixed hybrid finite element method and particle tracking algorithm.

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

The _cuDFNsys_ relies on several open packages: CUDA, Eigen, Gmsh, UMFPACK and HDF5.

Installations of these libraries could be disturbing, even you are familiar with Linux, CMake, C++ and so on. I cannot write a bash file that help you install all these libraries at one time. You should go to their homepages, and install them one by one. You may feel painful when install some libraries, for instance, the cuda package is difficult to install, and problem might happen occasionally. The Gmsh C++ API that _cuDFNsys_ relies on should support OCC mode, meaning that the occt library is required. Also, the installation of UMFPACK could be difficult. _Anyway, I am very willing to help you install them, but you should be familiar with Linux, CMake, C++ and so on_.

# Environment setup
After installation of above-mentioned packages, the _cuDFNsys_ should be linked to them. In directory 'Modules', .cmake files can be found, where the path of corresponding package should be added/set.

# Visualization

After simulation, _cuDFNsys_ outputs .h5 files, and you can visualize them by the generated .m file or .py file. 

Visulization with MATLAB is simple, you just need the license of MATLAB. 

With Python, you need to install the mayavi Python package. Sorry again for the inconvenience, the installation of mayavi could be difficult, but the mayavi is a very good visualization engine. The Matplotlib package is not good at 3D visualization, see https://matplotlib.org/2.2.2/mpl_toolkits/mplot3d/faq.html.

# Manual
No manual about the use of _cuDFNsys_ is available, because this is just for personal, scientific use right now. I may do that in the future.

# Directories

_CodingGuideline_: my coding guidline, e.g. how to name a member variable in a class

_Modules_: cmake script to find packages that the _cuDFNsys_ relies on in a Ubuntu system

_include_: header files containing declarations of functions/classes

_src_: source files containing definitions of functions/classes

_PercolationTest_, _TestFunctionsOnlyForDevlopment_, _TestResolutionEffect_ and other directories are all about examples (user's interfaces), with different purposes, meanwhile showing that how to call _cuDFNsys_ functions, with cmake and other tools.