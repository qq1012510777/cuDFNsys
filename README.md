# cuDFNsys

The _cuDFNsys_ is an open-source C++/CUDA library for DFN simulation, based on CUDA. It is also able to simulate static flow in DFN, based on mixed hybrid finite element method.

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

The _cuDFNsys_ relies on several packages: CUDA, Eigen, Gmsh, Matlab, UMFPACK, HDF5.

# Directories

_CodingGuideline_: my coding guidline, e.g. how to name a member variable in a class

_Modules_: cmake script to find packages that the _cuDFNsys_ relies on in a Ubuntu system

_include_: header files containing declarations of functions/classes

_src_: source files containing definitions of functions/classes

_PercolationTest_, _TestFunctionsOnlyForDevlopment_, _TestResolutionEffect_ and other directories are all about examples (user's interfaces), with different purposes, meanwhile showing that how to call _cuDFNsys_ functions, with cmake and other tools.