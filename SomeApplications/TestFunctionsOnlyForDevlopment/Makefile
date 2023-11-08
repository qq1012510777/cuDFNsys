ExeName=main
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
	$(NVCC) -DUSE_DOUBLES ./src/$(ExeName).cu $(cuDFNsysLibraryPath)/libcuDFNsys.a \
	-o ./$(ExeName) $(INCDIRS) $(Lib_) \
	-Xcudafe --display_error_number -Xcudafe --diag_suppress=3057 \
	-Xcudafe --diag_suppress=1301 \
	-Xcudafe --diag_suppress=3059 \
	-Xcudafe --diag_suppress=3060 \
	-Xcudafe --diag_suppress=3058 \
	-Xcudafe --diag_suppress=3056 \
	-arch=sm_60 -std=c++17 -rdc=true
clean:
	gio trash -f ./$(ExeName) *.h5 *.m *.py   log.txt