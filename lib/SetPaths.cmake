#***************************************************************************#
# cuDFNsys - simulating flow and transport in 3D fracture networks          #
# Copyright (C) 2022, Tingchang YIN, Sergio GALINDO-TORRES                  #
#                                                                           #
# This program is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU Affero General Public License as            #
# published by the Free Software Foundation, either version 3 of the        #
# License, or (at your option) any later version.                           #
#                                                                           #
# This program is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Affero General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU Affero General Public License  #
# along with this program.  If not, see <https://www.gnu.org/licenses/>.    #
#***************************************************************************#


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