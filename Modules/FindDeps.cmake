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

set (LIBs)
set (HEADERs)

# search FindXXX.cmake
# INCLUDE(${CUDFNSYS_ROOT}/Modules/FindMATLAB.cmake) # do not use MATLAB C++ API, because it is not free
INCLUDE(${CUDFNSYS_ROOT}/Modules/FindEIGEN.cmake)
INCLUDE(${CUDFNSYS_ROOT}/Modules/FindGMSH.cmake)
INCLUDE(${CUDFNSYS_ROOT}/Modules/FindUMFPACK.cmake)

FIND_PACKAGE (HDF5 COMPONENTS     HL)                       
FIND_PACKAGE (HDF5 COMPONENTS CXX HL)                       

OPTION(A_USE_OMP            "Use OpenMP  ?"                                        ON )

INCLUDE (FindOpenMP ) 
if(OPENMP_FOUND AND A_USE_OMP)
    ADD_DEFINITIONS (-DUSE_OMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    SET(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=${OpenMP_CXX_FLAGS}") 
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}") 
else(OPENMP_FOUND AND A_USE_OMP)
    if(A_USE_OMP)
        SET (MISSING "${MISSING} OpenMP")
    endif(A_USE_OMP)
endif(OPENMP_FOUND AND A_USE_OMP)

FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

if(HDF5_FOUND AND A_USE_HDF5)
    ADD_DEFINITIONS (-DH5_NO_DEPRECATED_SYMBOLS -DH5Gcreate_vers=2 -DH5Gopen_vers=2 -DUSE_HDF5)
	INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
	SET (LIBS ${LIBS} ${HDF5_LIBRARIES})
else(HDF5_FOUND AND A_USE_HDF5)
    if(A_USE_HDF5)
        SET (MISSING "${MISSING} HDF5")
    endif(A_USE_HDF5)
endif(HDF5_FOUND AND A_USE_HDF5)

# IF (MATLAB_FOUND)
#     SET (HEADERs ${HEADERs} ${MATLAB_INCLUDE_DIR})  
#     SET (LIBs ${LIBs} ${MATLAB_LIBRARIES}) 
# ELSE (MATLAB_FOUND)
#     SET (MISSING "${MISSING} MATLAB CXX API")   
# ENDIF (MATLAB_FOUND)

IF (EIGEN_FOUND)
    SET (HEADERs ${HEADERs} ${EIGEN_DENSE})  
ELSE (EIGEN_FOUND)
    SET (MISSING "${MISSING} EIGEN")   
ENDIF (EIGEN_FOUND)

IF (GMSH_FOUND)
    SET (HEADERs ${HEADERs} ${GMSH_INCLUDE_DIR})  
    SET (LIBs ${LIBs} ${GMSH_LIBRARIES}) 
ELSE (GMSH_FOUND)
    SET (MISSING "${MISSING} GMSH")   
ENDIF (GMSH_FOUND)

IF (UMFPACK_FOUND)
    SET (HEADERs ${HEADERs} ${UMFPACK_INCLUDE_DIRS})  
    SET (LIBs ${LIBs} ${UMFPACK_LIBRARIES}) 
ELSE (UMFPACK_FOUND)
    SET (MISSING "${MISSING} UMFPACK")   
ENDIF (UMFPACK_FOUND)

IF (HDF5_FOUND)
    SET (HEADERs ${HEADERs} ${HDF5_INCLUDE_DIR})  
    SET (LIBs ${LIBs} ${HDF5_LIBRARIES}) 
ELSE (HDF5_FOUND)
    SET (MISSING "${MISSING} HDF5")   
ENDIF (HDF5_FOUND)