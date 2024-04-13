# Manual for _cuDFNsys_

* Name: Tingchang YIN

* Institution: Westlake University, China

* Date: Oct. 19, 2023

* Update date: Apr. 13, 2024 

* Email: yintingchang@foxmail.com

**Now, _cuDFNsys_ has a simple GUI. Details are presented in Section 6.**

## 1. Compilation of a user interface code with cuDFNsys

Suppose that the _cuDFNsys_ library has been installed in `~/cuDFNsys/lib`, then `libcuDFNsys.a` can be seen in `~/cuDFNsys/lib`, which is a static library. How can one compile and install _cuDFNsys_ library (with nessesary dependencies, e.g., gmsh) can be seen in `README.md`.

Now, we can create a directory `CompilationTest` under `~/cuDFNsys`. We should test the compilation under this directory:
```
cd ~/cuDFNsys
mkdir CompilationTest
cd CompilationTest
touch Makefile
mkdir src
cd src
touch CompilationTest.cu
```

We can copy the following code to `~/cuDFNsys/CompilationTest/src/CompilationTest.cu`:
```
#include "cuDFNsys.cuh"
int main()
{
    int dev = 0;
    GPUErrCheck(cudaSetDevice(dev));
    cuDFNsys::Vector4<double> vec4 = cuDFNsys::MakeVector4(1.5, 1., 15., 0.);
    std::cout << "cuDFNsys::Vector4: " << vec4.x << ", " << vec4.y << ", " << vec4.z << ", " << vec4.w << ", " << std::endl;
    return 0;
}
```
The above code just generates a vector (type: `cuDFNsys::Vector4<double>`) named `vec4` of four elements (double precision), which can be accessed by `vec4.x`, `vec4.y`, `vec4.z` and `vec4.w`. This vector is actually `double3` in cuda!!!

To compile this code linking to `libcuDFNsys.a`, we can copy the following `Make` script to `~/cuDFNsys/CompilationTest/Makefile`:
```
ExeName=CompilationTest
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
```

Note that, if necessary, some variables of this script should be changed for your Ubuntu, namely,
```
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
```
The above variables are the paths to different tools, CXX head files, and CXX libraries that _cuDFNsys_ depends. For example, you can change the variable `GmshIncludePath` like
```
GmshIncludePath=$(HOME)/gmsh/include
```
if you installed `gmsh-dev` under `~/gmsh` in your Ubuntu.

Now, we can comile `CompilationTest.cu` by:
```
cd ~/cuDFNsys/CompilationTest
make
```
If no error happens, then an executable file `main` should appear in `~/cuDFNsys/CompilationTest`. We just run it by
```
./main
```
and we should see
```
cuDFNsys::Vector4: 1.5, 1, 15, 0,
```
printed on the terminal window.

If errors happen when you run `./main`, e.g., `error while loading shared libraries: libgmsh.so: cannot open shared object file: No such file or directory`, two enviromental variables should be added to `~/.bashrc` as follows:
```
export LD_LIBRARY_PATH=path-to-gmsh-library:$LD_LIBRARY_PATH
export LIBRARY_PATH=path-to-gmsh-library:$LIBRARY_PATH
# change path-to-gmsh-library to the dynamic gmsh library in the computer
```
then update `~/.bashrc` by
```
source ~/.bashrc
```
This is a runtime error, that the computer does not know where the library `gmsh` is when running `./main`.

## 2. Use an executable file
To users who do not want to write a `.cpp` file, an executable file - `cuDFNsys_exe`, can be used.

First, the executable should be compiled.
```
cd ~/cuDFNsys/bin
make
```

Then, you can see two files, `cuDFNsys_exe` and `Transform2DH5ParticleDataTo3D`, generated in the current directory.

`cuDFNsys_exe` needs two input `.csv` files to get fracture attributes and other parameters.

To generate DFNs, I provide two `.csv` files.

### 2.1 csv to generate a deterministic DFN
```
#InputParametersForDeterministicDFN.csv
IfStochastic,0,
DomainSizeX,30,
DomainDimensionRatio,1,1,1,,,,,
Percolation_direction,2,,,,,,,
NumFractures,4,
Fracture_1, 1,0,0,0,0,0,10,0.2,5e-4,0,1,1,
Fracture_2, 0,0,1,0,0,0,55,0.1,7e-9,,,,
Fracture_3, 0,1,0,0,0,0,55,0.3,6e-7,
Fracture_4, 1,1,1,10,10,-25,4,0.3,6e-7,

```
This one means a deterministic DFN, defined by users.
`IfStochastic,0,` denotes that this DFN is deterministic.
`DomainSizeX,30,` denotes that the length of the side in the `x` direction is 30 m.
`DomainDimensionRatio,1,1,1,,,,,` denotes that the dimension ratio of the DFN domain is 1 : 1 : 1, thus the size of the domain is 30m, 30m, 30m
`Percolation_direction,2,,,,,,,` means that the program will check if the DFN is percolative in `z` direction
`NumFractures,4,` means that there are four fractures
`Fracture_1, 1,0,0,0,0,0,10,0.2,5e-4,0,1,1,`: there are 12 numbers separated by commas, representing the fracture attribues of a fracture. The first three `1,0,0` means the normal vector of the fracture plane. Then, the next three `0,0,0` is the center coodinate. The seventh `10` is the radius $R$ of the fracture. The eighth `0.2` is the $\beta$ value of fractures, the ninth `5e-4` is the $\gamma$ value. In _cuDFNsys_, the conductivity $k_f = b_f^3 / 12 =  \left[(\gamma R ^\beta) \right] ^ 3 / 12$, where $b_f$ is the aperture of fractures. The last three `0,1,1` denotes the coordinates of one vertices of the square fracture. If the last three number is `0,0,0,` or `,,,` (no numbers), then the program will randomly determine a vertex in the fracture plane with $R$.

### 2.2 csv to generate a random DFN
```
##InputParametersForStochasticDFN.csv
IfStochastic,1,
DomainSizeX,30,
DomainDimensionRatio,1,1,2,,,,,
Percolation_direction,2,,,,,,,
NumFractureGroups,3,,,,,,,
NumFractureEachGroup,70,50,30,,,,,
KappaValues,20,20,25,,,,,
MeanOrientationOfFisherDistribution,0,0,1,1,0,0,0,1,0,
ModeOfSizeDistribution,0,1,2,,,,,
SizeDistributionParameters,1.5,1,15,0,8.5,5.5,1,15,1,15,0,0,
Beta,0.2,0.3,0.4,
Gamma,2.0e-5,3.0e-6,4e-5,
```

`IfStochastic,1,` denotes a stochastic DFN.
`NumFractureGroups,3,,,,,,,` denotes that there three fracture groups.
`KappaValues,20,20,25,,,,,`: the $\kappa$ values (Fisher constants) of the three groups.
`MeanOrientationOfFisherDistribution,0,0,1,1,0,0,0,1,0,`: mean orientations of the three groups in vector forms.
`ModeOfSizeDistribution,0,1,2,,,,,`: the size distribution patterns of each fracture Group. When the mode is `0`, the size distribution is a power-law. When the mode  is `1`, the size distribution is lognormal. When the mode  is `2`, the size distribution is uniform. When the mode  is `3`, the size distribution is mono-sized.
`SizeDistributionParameters`: the statistic parameters of the size distributions of each group. Four numbers is a set of statistic parameters for a size distribution of a group. `1.5,1,15,0` means, as the distribution of first group is a power-law, the power-law exponent is 1.5, minimum $R$ is 1 m, maximum is 15 m. As the second group has a longnormal size distribution, `8.5,5.5,1,15` means that the mean is 8.5, variance is 5.5, min R is 1, max $R$ is 15 m. The this group has a uniform size distribution, thus `1,15,0,0,` means that min $R$ is 1, max $R$ is 15 m.
`Beta,0.2,0.3,0.4,` and `Gamma,2.0e-5,3.0e-6,4e-5,` mean the $\beta$ and $\gamma$ values of each group, respectively.

### 2.3. Run a DFN generation and mesh and flow simulation
One can run the program by `./cuDFNsys_exe`. Then the terminal shows
```
Creating a new DFN. Please provide the name of the .csv file, without the suffix .csv!
```
One can input the name of `.csv`, either `InputParametersForDeterministicDFN` or `InputParametersForStochasticDFN`.

Then, the program will load parameters for mesh and flow simulation, in `Mesh_parameters.csv` and `PT_parameters.csv`, respectively, e.g., the minimum and maximum expected finite element sizes, the inlet and outlet head values, and so on.

### 2.4. Particle tracking
```
### PT_parameters.csv
NumParticles,20000,
NumTimeSteps,2000,
MolecularDiffusion,1.39715e-12,
DeltaT, 76491223.69119304,
InjectionMethod_initialcondition,Flux-weighted,
If_OutputAllPTInformationOrFPTCurve, 1,
SpacingOfControlPlanes, 30,
IfOutputVarianceOfDisplacementsEachStep, 1,
IfInjectAtCustomedPlane, 1,
CustomedPlaneInjection, 23,
IfUseFluxWeightedOrEqualProbableMixingIntersection, 1,
```
The parameters of particle tracking are defined in a `.csv`.

`NumParticles` is the number of particles injected into DFNs.
`InjectionMethod_initialcondition` is the initial condition, which can be `Flux-weighted`, `Resident`, or `Point`.
`If_OutputAllPTInformationOrFPTCurve` can be 0 or 1. The value `1` means particle coordinates of all steps are outputted.
`SpacingOfControlPlanes` is the spacing of control planes.
`IfOutputVarianceOfDisplacementsEachStep` can be 0 or 1. The value `1` means that the variance of particle displacements is outputted for each step.
`IfInjectAtCustomedPlane` can be 0 or 1. The value `1` means that the injection plane is user-defined.
`CustomedPlaneInjection`: if `IfInjectAtCustomedPlane == 1`, then `IfInjectAtCustomedPlane` is the location of injection plane.
`IfUseFluxWeightedOrEqualProbableMixingIntersection`: `0` is equiproble mixing rule at intersections. `1`: outgoing-flux-weighted mixing rule at intersections.

### 2.5. Visualization
After the simulation, one can see some `.m` or `.py` files in the current directory. One can run them by MATLAB or python to visualize DFN, DFN mesh and DFN hydraulic head field.

### 2.6. Data of DFN, mesh, and flow
These will be discussed in the next section.

## 3. Data files
By running `cuDFNsys_exe`, three `.h5` files are outputted, e.g., `Class_DFN.h5`, `Class_MESH.h5` and `Class_FLOW.h5`.

One can open them to see DFN parameters.

The first one `Class_DFN.h5` has the following structure:
```
/                        Group
/Cluster_1               Dataset {338}
/ *******
/Cluster_9               Dataset {1}

/DomainDimensionRatio    Dataset {3}
/Fracture_1              Group
/Fracture_1/Center       Dataset {1, 3}
/Fracture_1/Conductivity Dataset {1}
/Fracture_1/ConnectModelSurf Dataset {1, 6}
/Fracture_1/NormalVec    Dataset {1, 3}
/Fracture_1/NumVertsTruncated Dataset {1}
/Fracture_1/Radius       Dataset {1}
/Fracture_1/Verts3D      Dataset {3, 4}
/Fracture_1/Verts3DTruncated Dataset {3, 8}

/IfPeriodic              Dataset {1}
/Intersections           Dataset {8, 1656}
/L                       Dataset {1}
/NumClusters             Dataset {1}
/NumFractures            Dataset {1}
/PercoDir                Dataset {1}
/PercolationClusters     Dataset {1}
```
where `Dataset {*}` is the dimension of the dataset. One can what each dataset is thourgh their names. Note that for other DFNs, the dimension of some datasets will be different. 

`Cluster_1` is the ID of fractures belonging to the first cluster.

Note that `./Fracture_1/NumVertsTruncated` is the number of vertices after the truncation of fractures. `Verts3DTruncated` is the vertices' coordinates after the truncation.

The second one `Class_MESH.h5` has the following structure:
```
/Fracs_percol            Dataset {86}
/InletTraceLength        Dataset {1}
/MeanGridSize            Dataset {1}
/NumElements             Dataset {1}
/NumFractures            Dataset {1}
/NumInletEdges           Dataset {1}
/NumNodes                Dataset {1}
/NumOutletEdges          Dataset {1}
/OutletTraceLength       Dataset {1}
/group_mesh              Group
/group_mesh/Coordinate2D Dataset {6, 44111}
/group_mesh/Coordinate3D Dataset {3, 22879}
/group_mesh/Dir          Dataset {1}
/group_mesh/EdgeAttri    Dataset {6, 44111}
/group_mesh/Element2D    Dataset {3, 44111}
/group_mesh/Element3D    Dataset {3, 44111}
/group_mesh/ElementFracTag Dataset {44111}
/group_mesh/FracID       Dataset {86}
/group_mesh/InletEdgeNOLen Dataset {2, 163}
/group_mesh/MeshSuccess  Dataset {1}
/group_mesh/NumElementFrac Dataset {86}
/group_mesh/NumInteriorEdges Dataset {1}
/group_mesh/NumNeumannEdges Dataset {1}
/group_mesh/OutletEdgeNOLen Dataset {2, 143}
```
where `Fracs_percol` is the tag/ID of fractures that belong to the percolation cluster. `InletTraceLength` is the total length of traces at the inlet plane. `MeanGridSize` is the mean area of finite elements. `NumElements` is the number of finite elements. `NumInletEdges` is the number of traces on inlet plane. `NumNodes` is the number of nodes of the unstructured mesh. `/group_mesh/Coordinate3D` is the 3D coordinates of all nodes. `Dir` is the percolation direction. `Element3D` is the structure of finite elements, showing us how a finite element consists of which nodes. `ElementFracTag` is the fracture ID for each finite element.

The third one `Class_FLOW.h5` has the following structure:
```
/Dir                     Dataset {1}
/InletLength             Dataset {1}
/InletP                  Dataset {1}
/MaxVelocity             Dataset {1}
/MeanVelocity            Dataset {1}
/MuOverRhoG              Dataset {1}
/OutletLength            Dataset {1}
/OutletP                 Dataset {1}
/Permeability            Dataset {1}
/PressureEles            Dataset {1, 44111}
/PressureInteriorEdge    Dataset {1, 66765}
/QError                  Dataset {1}
/QIn                     Dataset {1}
/QOut                    Dataset {1}
/TripletTime             Dataset {1}
/VelocityNormalScalarSepEdges Dataset {1, 132333}
```
where `MaxVelocity` is the maximum flow rate $\left[LT^{-1}\right]$, `MuOverRhoG` is the viscosity $\mu$ over the product of fluid density $\rho_g$ and gravitational acceleration $g$, `Permeability` is the Darcian permeability, `PressureEles` is the hydraulic head of each element, `VelocityNormalScalarSepEdges` is the volumatric flux rate per width of each edge (each element has three edges).

Finally, there are some `.h5` files storing particle tracking informations. `cuDFNsys/bin/Dispersion_MeanSquareDisplacement.h5` stores the variances of particle displacements in the $x$, $y$ and $z$ directions at each time step. `cuDFNsys/bin/ParticlePositionResult/ParticlePositionInit.h5` stores the initial position of particles' coodinates (which is in the 2D format). `cuDFNsys/bin/ParticlePositionResult/ParticlePositionBlock000000000_n_.h5` (`_n_` could be a number, e.g., 1, 2, 3 or 4) stores the positions of particles at each step. Note that the 'cuDFNsys_exe' program will ask the user that if they want to transform 2D coodinates to 3D coodinates. After the transformation, the dispersion can be visualized.

## 4. Quickstart examples: write a main.cpp
The whole code is shown in `cuDFNsys/QuickStartGuide/src/QuickStartGuide.cpp`. To compile this, `cd cuDFNsys/QuickStartGuide` and `make`, then the executable program `QuickStartGuide` is generated in the current directory.

First, if one want to generate a stochastic DFN, the `time` function is required.
```
time_t t;
time(&t);
```

Next, one can create four empty C++ classes to represent DFN, mesh, flow and particle tracking:
```
cuDFNsys::DFN<double> my_dfn;
cuDFNsys::MeshDFN<double> meshGen;
cuDFNsys::FlowDFN<double> flowDFN;
cuDFNsys::PTDFN<double> particleTracking;
```

To generate a stochastic DFN, one can load input parameters from a `.csv` file:
```
my_dfn.RandomSeed = (unsigned long)t;
my_dfn.LoadDFNFromCSV("InputParametersForStochasticDFN");
```
Note that `RandomSeed` aims at stochastic DFNs at different time. Right now, a DFN is created. The intersections and clusters (including spanning clusters) should be identified, which can be done by
```
my_dfn.IdentifyIntersectionsClusters(true); // true here means consideration of truncated fractures
```
This DFN now can be visualized. The visulization can be performed by python or MATLAB. One can output the visualization files  by
```
my_dfn.Visualization("DFN_VISUAL", "DFN_VISUAL", "DFN_VISUAL", true, true,
                         true, true);
```
The above function generates `.m`, `.py` and `.h5` files. For instance, one can visualize the DFN by running `DFN_VISUAL.m` in MATLAB.

Based on the DFN, mesh, flow simulation and particle tracking can be done by:
```
meshGen.LoadParametersFromCSV("Mesh_parameters");
meshGen.MeshGeneration(my_dfn);
meshGen.Visualization(my_dfn, "DFN_MESH_VISUAL", "DFN_MESH_VISUAL",
                          "DFN_MESH_VISUAL", true, true);

flowDFN.LoadParametersFromCSV("Flow_parameters");
flowDFN.FlowSimulation(my_dfn, meshGen);
flowDFN.Visualization(my_dfn, meshGen, "DFN_FLOW_VISUAL", "DFN_FLOW_VISUAL",
                          "DFN_FLOW_VISUAL");		

particleTracking.LoadParametersFromCSV("PT_parameters");
particleTracking.ParticleTracking(my_dfn, meshGen, flowDFN);
particleTracking.Visualization(my_dfn, meshGen, flowDFN,
                                   "DFN_DISPERSION_VISUAL",
                                   "DFN_DISPERSION_VISUAL", "DFN_FLOW_VISUAL");
```
All input parameters are shown in the corresponding `.csv` files.

Finally, the DFN, mesh and flow data can be stored (for re-run particle tracking or other purpose):
```
my_dfn.StoreInH5("Class_DFN");
meshGen.StoreInH5("Class_MESH");
flowDFN.StoreInH5("Class_FLOW");
```
The three h5 files can be loaded again for run more particle tracking steps (details are shown in the next section).

To show the animation of particle tracking, one can run `./Transform2DH5ParticleDataTo3D 0 DFN_MESH_VISUAL.h5` in the terminal. Then run `DFN_DISPERSION_VISUAL.m` in MATLAB.

## 5. How to re-run a DFN
In the last section, the generation of a DFN, mesh, flow simulation, and particle tracking have been presented. A question arises here: how can one run more particle tracking steps, based on the existing results?

A re-run code can be written. The whole code can be seen in `cuDFNsys/QuickStartGuide/src/ReRunPT.cpp`.

First, again, the four classes should be created:
```
    cuDFNsys::DFN<double> my_dfn;
    cuDFNsys::MeshDFN<double> meshGen;
    cuDFNsys::FlowDFN<double> flowDFN;
    cuDFNsys::PTDFN<double> particleTracking;
```
Then, the previous data files can be loaded by 
```
    my_dfn.LoadClassFromH5("Class_DFN");
    meshGen.LoadClassFromH5("Class_MESH");
    flowDFN.LoadClassFromH5("Class_FLOW");
```
Finally, one just need to load the particle tracking parameters again, in which only the number of times step is taken, because other parameters, like diffusion coefficients, have been fixed.
```
    particleTracking.LoadParametersFromCSV("PT_parameters");
    particleTracking.ParticleTracking(my_dfn, meshGen, flowDFN);
```

## 6. GUI
Now, _cuDFNsys_ has a simple GUI, which can be run on Ubuntu or Windows WSL. To use this, make sure the _cuDFNsys_ library has been compiled: the file `libcuDFNsys.a` should appear in `~/cuDFNsys/lib`.

Before compiling the GUI, the `Path.mk` under `~/cuDFNsys` should be changed, which is a included file of Makefile, containing paths of different denpendencies. These paths should be changed if you installed them at different paths. Details can be found in Section 1.

First, `cd cuDFNsys/GUI`, and then `make`. The `make` command here will compile four executables in the current directory: `DFN_Gen`, `DFN_Mesh`, `DFN_Flow` and `Transform2DH5ParticleDataTo3D`. So, once `make` command finished, you can see the four executables in the `cuDFNsys/GUI` directory. How to build _cuDFNsys_ can be found in README.md.

Now, we can run the GUI by `python3 src/cuDFNsysGUI.py`, as shown below.

cuDFNsys GUI can be run under other directories. This can be done by adding one line to the end of `~/.bashrc`:
```
echo "alias cuDFNsys='python3 ~/cuDFNsys/GUI/src/cuDFNsysGUI.py'" >> ~/.bashrc
``` 
One can run `source ~/.bashrc` to update it.  Then one can just input `cuDFNsys` to run GUI under any directory.

<p align="center">
  <img width="500" src="https://github.com/qq1012510777/cuDFNsys/blob/main/Manual/GUI.png">
</p>
<p align="center">
    <em>Figure 1. The GUI of cuDFNsys. </em>
</p>

Before running a DFN, a project (i.e., a directory) should be established first, by clicking `File/New project`. One can also open an existing project by `Open project`. By the way, the root path of _cuDFNsys_ is defaultly set to be under home directory. If one want to change it, just click `Set cuDFNsys root`.

To generate a DFN, one can click `DFN generation`. One can generate a Stochastic DFN, or a deterministic DFN, or a DFN contaning both stochastic and deterministic fractures. To visualize the DFN, click `DFN generation/Use matplotlib` or `DFN generation/Use mayavi`.

The mesh, flow and particle tracking simulations are also similar. When click corresponding buttons, one should fill parameters to entries, the meaning of the entries can be obtained by click red `Help` buttons.

The GUI also supports Monte Carlo iterations under the same condition (stochastic DFNs), just for DFN and mesh generation and flow simulation.

For visualizations, the GUI will use Matplotlib or Mayavi (see README.md for its installation).