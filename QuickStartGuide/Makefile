ExeName=QuickStartGuide
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

INCDIRS=-I ../include \
		-I $(Hdf5IncludePath) \
		-I $(GmshIncludePath) \
		-I $(EigenIncludePath) \
		-I /usr/local/include \
		-I $(UmfpackIncludePath)
Lib_= -Xcompiler=-fopenmp \
	-lpthread \
	-lm \
	-lcudadevrt \
	-lcudart_static \
	-lrt \
	-lpthread \
	-ldl \
	-lgmsh \
	-lumfpack \
	-lamd \
	-lhdf5_cpp \
	-lhdf5 \
	-lsz \
	-lz \
	-ldl \
	-L$(GmshLibraryPath) \
	-L$(UmfpackLibraryPath) \
	-L$(Hdf5LibraryPath)
# -L/usr/lib/x86_64-linux-gnu 
# -L"/usr/local/cuda-11.0/targets/x86_64-linux/lib/stubs" \
# -L"/usr/local/cuda-11.0/targets/x86_64-linux/lib"

# -l$(HOME)/pkg/gmsh-4.8.4-source/MY_GMSH/lib/libgmsh.so \
# -l$(HOME)/pkg/SuiteSparse-master/lib/libumfpack.so \
# -l$(HOME)/pkg/SuiteSparse-master/lib/libamd.so \
# -l/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a \
# -l/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so \
# -l/usr/lib/x86_64-linux-gnu/libsz.so \
# -l/usr/lib/x86_64-linux-gnu/libz.so \
# -l/usr/lib/x86_64-linux-gnu/libdl.so  

all: 
	$(NVCC) -DUSE_DOUBLES ./src/$(ExeName).cu ../lib/libcuDFNsys.a \
	-o ./$(ExeName) $(INCDIRS) $(Lib_) \
	-Xcudafe --display_error_number -Xcudafe --diag_suppress=3057 \
	-Xcudafe --diag_suppress=1301 \
	-Xcudafe --diag_suppress=3059 \
	-Xcudafe --diag_suppress=3060 \
	-Xcudafe --diag_suppress=3058 \
	-arch=sm_60 -std=c++17 -rdc=true

	$(NVCC) -DUSE_DOUBLES ../Transform2DH5ParticleDataTo3D/src/Transform2DH5ParticleDataTo3D.cu ../lib/libcuDFNsys.a \
	-o ./Transform2DH5ParticleDataTo3D $(INCDIRS) $(Lib_) \
	-Xcudafe --display_error_number -Xcudafe --diag_suppress=3057 \
	-Xcudafe --diag_suppress=1301 \
	-Xcudafe --diag_suppress=3059 \
	-Xcudafe --diag_suppress=3060 \
	-Xcudafe --diag_suppress=3058 \
	-arch=sm_60 -std=c++17 -rdc=true
clean: