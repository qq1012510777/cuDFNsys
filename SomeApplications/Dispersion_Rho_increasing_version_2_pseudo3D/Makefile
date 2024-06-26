
include ../Path.mk

INCDIRS=-I $(cuDFNsysIncludePath)\
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

all: main DFN_Gen DFN_Mesh DFN_Flow DFN_PT

.PHONY: all main DFN_Gen DFN_Mesh DFN_Flow DFN_PT

main: ./src/main.cu
	$(NVCC) -DUSE_DOUBLES ./src/main.cu $(cuDFNsysLibraryPath)/libcuDFNsys.a \
	-o ./main $(INCDIRS) $(Lib_) \
	-Xcudafe --display_error_number -Xcudafe --diag_suppress=3057 \
	-Xcudafe --diag_suppress=1301 \
	-Xcudafe --diag_suppress=3059 \
	-Xcudafe --diag_suppress=3060 \
	-Xcudafe --diag_suppress=3058 \
	-Xcudafe --diag_suppress=3056 \
	-arch=sm_60 -std=c++17 -rdc=true


DFN_Gen: ./src/DFN_Gen.cu 
	$(NVCC) -DUSE_DOUBLES ./src/DFN_Gen.cu $(cuDFNsysLibraryPath)/libcuDFNsys.a \
	-o ./DFN_Gen $(INCDIRS) $(Lib_) \
	-Xcudafe --display_error_number -Xcudafe --diag_suppress=3057 \
	-Xcudafe --diag_suppress=1301 \
	-Xcudafe --diag_suppress=3059 \
	-Xcudafe --diag_suppress=3060 \
	-Xcudafe --diag_suppress=3058 \
	-arch=sm_60 -std=c++17 -rdc=true


DFN_Mesh: ./src/DFN_Mesh.cu
	$(NVCC) -DUSE_DOUBLES ./src/DFN_Mesh.cu $(cuDFNsysLibraryPath)/libcuDFNsys.a \
	-o ./DFN_Mesh $(INCDIRS) $(Lib_) \
	-Xcudafe --display_error_number -Xcudafe --diag_suppress=3057 \
	-Xcudafe --diag_suppress=1301 \
	-Xcudafe --diag_suppress=3059 \
	-Xcudafe --diag_suppress=3060 \
	-Xcudafe --diag_suppress=3058 \
	-arch=sm_60 -std=c++17 -rdc=true

DFN_Flow: ./src/DFN_Flow.cu
	$(NVCC) -DUSE_DOUBLES ./src/DFN_Flow.cu $(cuDFNsysLibraryPath)/libcuDFNsys.a \
	-o ./DFN_Flow $(INCDIRS) $(Lib_) \
	-Xcudafe --display_error_number -Xcudafe --diag_suppress=3057 \
	-Xcudafe --diag_suppress=1301 \
	-Xcudafe --diag_suppress=3059 \
	-Xcudafe --diag_suppress=3060 \
	-Xcudafe --diag_suppress=3058 \
	-arch=sm_60 -std=c++17 -rdc=true

DFN_PT: ./src/DFN_PT.cu
	$(NVCC) -DUSE_DOUBLES ./src/DFN_PT.cu $(cuDFNsysLibraryPath)/libcuDFNsys.a \
	-o ./DFN_PT $(INCDIRS) $(Lib_) \
	-Xcudafe --display_error_number -Xcudafe --diag_suppress=3057 \
	-Xcudafe --diag_suppress=1301 \
	-Xcudafe --diag_suppress=3059 \
	-Xcudafe --diag_suppress=3060 \
	-Xcudafe --diag_suppress=3058 \
	-arch=sm_60 -std=c++17 -rdc=true

clean: