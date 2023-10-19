ExeName=CompilationTest
NVCC=/usr/local/cuda/bin/nvcc
# include paths for headers
cuDFNsysIncludePath=$(HOME)/cuDFNsys/include
Hdf5IncludePath=/usr/lib/x86_64-linux-gnu/hdf5/serial/include
GmshIncludePath=$(HOME)/pkg/gmsh-4.8.4-source/MY_GMSH/include
EigenIncludePath=$(HOME)/pkg/eigen
UmfpackIncludePath=$(HOME)/pkg/SuiteSparse-master/include
# library paths
cuDFNsysLibraryPath=$(HOME)/cuDFNsys/lib
GmshLibraryPath=$(HOME)/pkg/gmsh-4.8.4-source/MY_GMSH/lib
UmfpackLibraryPath=$(HOME)/pkg/SuiteSparse-master/lib
Hdf5LibraryPath=/usr/lib/x86_64-linux-gnu/hdf5/serial

INCDIRS=-I $(cuDFNsysIncludePath) \
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
all: 
	$(NVCC) -DUSE_DOUBLES ./src/$(ExeName).cu $(cuDFNsysLibraryPath)/libcuDFNsys.a \
	-o ./main $(INCDIRS) $(Lib_) \
	-Xcudafe --display_error_number -Xcudafe --diag_suppress=3057 \
	-Xcudafe --diag_suppress=1301 \
	-Xcudafe --diag_suppress=3059 \
	-Xcudafe --diag_suppress=3060 \
	-Xcudafe --diag_suppress=3058 \
	-arch=sm_60 -std=c++17 -rdc=true
clean: