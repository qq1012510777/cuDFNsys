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

cmake_minimum_required (VERSION 3.5)

project (benchmark)

enable_language(CUDA)

set (CUDFNSYS_ROOT $ENV{HOME}/cuDFNsys)

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin) ## bin is the EXECUTABLE_OUTPUT_PATH

set (SRC_LIST ${PROJECT_SOURCE_DIR}/src/main.cu)

include_directories(${CUDFNSYS_ROOT}/include)

INCLUDE(${CUDFNSYS_ROOT}/Modules/FindDeps.cmake)

include_directories(${HEADERs})

ADD_DEFINITIONS(-std=c++17)  
ADD_DEFINITIONS(-DUSE_DOUBLES) ## use double data type

set (SRC_LIST ${PROJECT_SOURCE_DIR}/src/main.cu)

file(GLOB_RECURSE SRC_CU_PATH ${CUDFNSYS_ROOT}/src/*.cu)

ADD_EXECUTABLE(main
    ${SRC_LIST} ${SRC_CU_PATH}
)

TARGET_LINK_LIBRARIES( main ${CUDA_LIBRARIES})
TARGET_LINK_LIBRARIES( main ${LIBs})

SET (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe --display_error_number")
SET (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe --diag_suppress=3057")
SET (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe --diag_suppress=1301")
SET (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe --diag_suppress=3059")
SET (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_60 -std=c++17 -rdc=true -O2 -DNDEBUG")