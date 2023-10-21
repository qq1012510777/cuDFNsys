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

cmake_minimum_required (VERSION 3.10)

project (generate_lib)

INCLUDE(${PROJECT_SOURCE_DIR}/SetPaths.cmake)

# find_package(CUDA)
# if (NOT CUDA_FOUND)
#     message(STATUS "CUDA not found. Project will not be built.")
# endif(NOT CUDA_FOUND)
if(NOT CMAKE_CUDA_COMPILER)
    set (CMAKE_CUDA_COMPILER ${NVCC_PATH})
ENDIF()
enable_language(CUDA)

set (CUDFNSYS_ROOT $ENV{HOME}/cuDFNsys)

include_directories(${CUDFNSYS_ROOT}/include)

INCLUDE(${CUDFNSYS_ROOT}/Modules/FindDeps.cmake)

include_directories(${HEADERs})

# ADD_DEFINITIONS(-std=c++17)  

file(GLOB_RECURSE SRC_CU_PATH ${CUDFNSYS_ROOT}/src/*.cu)

# SET (CMAKE_VERBOSE_MAKEFILE TRUE)
SET (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe --display_error_number")
SET (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe --diag_suppress=3057")
SET (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe --diag_suppress=1301")
SET (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe --diag_suppress=3059")
SET (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe --diag_suppress=20014")
SET (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_60 -std=c++17 -rdc=true")

# SET (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcudafe --display_error_number")
# SET (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcudafe --diag_suppress=3057")
# SET (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcudafe --diag_suppress=1301")
# SET (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcudafe --diag_suppress=3059")
# SET (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_60 -std=c++17 -rdc=true -O3 -Xcompiler -fPIC")

#add_library (cuDFNsys_shared SHARED ${SRC_CU_PATH})
add_library(cuDFNsys_static STATIC ${SRC_CU_PATH})

#set_target_properties(cuDFNsys_shared PROPERTIES OUTPUT_NAME "cuDFNsys")
set_target_properties(cuDFNsys_static PROPERTIES OUTPUT_NAME "cuDFNsys")

set (LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

# TARGET_LINK_LIBRARIES( cuDFNsys_shared ${CUDA_LIBRARIES} ${LIBs})

TARGET_LINK_LIBRARIES( cuDFNsys_static ${CUDA_LIBRARIES} ${LIBs})

# set_target_properties(cuDFNsys_shared
#                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

set_target_properties(cuDFNsys_static
                      PROPERTIES CUDA_SEPARABLE_COMPILATION ON)