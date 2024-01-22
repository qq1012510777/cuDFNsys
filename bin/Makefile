ExeName=cuDFNsys_exe
# NVCC path
NVCC=/usr/lib/nvidia-cuda-toolkit/bin/nvcc
# include paths for headers
cuDFNsysIncludePath=$(HOME)/cuDFNsys/include
Hdf5IncludePath=/usr/include/hdf5/serial
GmshIncludePath=usr/include
EigenIncludePath=/usr/include
UmfpackIncludePath=/usr/include/suitesparse
# library paths
cuDFNsysLibraryPath=$(HOME)/cuDFNsys/lib
GmshLibraryPath=/usr/lib/x86_64-linux-gnu
UmfpackLibraryPath=/usr/lib/x86_64-linux-gnu
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

Transform2DH5ParticleDataTo3Dname=Transform2DH5ParticleDataTo3D

all:  
	$(NVCC) -DUSE_DOUBLES ./Exe_src/main.cu $(cuDFNsysLibraryPath)/libcuDFNsys.a \
	-o ./$(ExeName) $(INCDIRS) $(Lib_) \
	-Xcudafe --display_error_number -Xcudafe --diag_suppress=3057 \
	-Xcudafe --diag_suppress=1301 \
	-Xcudafe --diag_suppress=3059 \
	-Xcudafe --diag_suppress=3060 \
	-Xcudafe --diag_suppress=3058 \
	-arch=sm_60 -std=c++17 -rdc=true


	$(NVCC) -DUSE_DOUBLES ../SomeApplications/Transform2DH5ParticleDataTo3D/src/Transform2DH5ParticleDataTo3D.cu $(cuDFNsysLibraryPath)/libcuDFNsys.a \
	-o ./Transform2DH5ParticleDataTo3D $(INCDIRS) $(Lib_) \
	-Xcudafe --display_error_number -Xcudafe --diag_suppress=3057 \
	-Xcudafe --diag_suppress=1301 \
	-Xcudafe --diag_suppress=3059 \
	-Xcudafe --diag_suppress=3060 \
	-Xcudafe --diag_suppress=3058 \
	-arch=sm_60 -std=c++17 -rdc=true

clean:
	gio trash -f ./*.h5 ./*.py ./*.m ./ParticlePositionResult 