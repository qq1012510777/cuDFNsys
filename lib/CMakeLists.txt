cmake_minimum_required (VERSION 3.7)

project (generate_lib)

# find_package(CUDA)
# if (NOT CUDA_FOUND)
#     message(STATUS "CUDA not found. Project will not be built.")
# endif(NOT CUDA_FOUND)

enable_language(CUDA)

set (CUDFNSYS_ROOT $ENV{HOME}/Desktop/cuDFNsys_II)

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
SET (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_60 -std=c++17 -rdc=true")

# SET (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcudafe --display_error_number")
# SET (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcudafe --diag_suppress=3057")
# SET (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcudafe --diag_suppress=1301")
# SET (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcudafe --diag_suppress=3059")
# SET (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_60 -std=c++17 -rdc=true -O3 -Xcompiler -fPIC")

add_library (cuDFNsys_shared SHARED ${SRC_CU_PATH})
add_library(cuDFNsys_static STATIC ${SRC_CU_PATH})

set_target_properties(cuDFNsys_shared PROPERTIES OUTPUT_NAME "cuDFNsys")
set_target_properties(cuDFNsys_static PROPERTIES OUTPUT_NAME "cuDFNsys")

set (LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

TARGET_LINK_LIBRARIES( cuDFNsys_shared ${CUDA_LIBRARIES} ${LIBs})

TARGET_LINK_LIBRARIES( cuDFNsys_static ${CUDA_LIBRARIES} ${LIBs})