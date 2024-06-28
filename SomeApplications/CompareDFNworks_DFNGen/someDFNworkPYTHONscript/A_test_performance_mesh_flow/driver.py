from pydfnworks import *
import os
import numpy as np
import h5py
import time

# import random

src_path = os.getcwd()

dfnFlow_file = src_path + "/dfn_explicit.in"
dfnTrans_file = src_path + "/PTDFN_control.dat"

n_cpu = 10
n_MC = 30

DFNGen_Time = np.zeros((n_cpu, n_MC))
DFNMesh_Time = np.zeros((n_cpu, n_MC))
DFNFlow_Time = np.zeros((n_cpu, n_MC))

DFN_mesh_flow_total_time = np.zeros((n_cpu, 1))

Index_iiii = 0
for i in range(n_cpu):

    start_time_iiiii = time.time()

    for j in range(n_MC):
        jobname = src_path + "/output_CPUThreads_" + str(i + 1) + "_MC_NO_" + str(j + 1)
        DFN = DFNWORKS(
            jobname, dfnFlow_file=dfnFlow_file, dfnTrans_file=dfnTrans_file, ncpu=i + 1
        )
        # sdsf_seed = int(time.time()) ^ int.from_bytes(os.urandom(4), byteorder='big')
        DFN.params["domainSize"]["value"] = [30, 30, 30]
        DFN.params["domainSizeIncrease"]["value"] = [0, 0, 0]
        DFN.params["h"]["value"] = 0.1
        DFN.params["stopCondition"]["value"] = 0
        DFN.params["nPoly"]["value"] = 90
        DFN.params["seed"]["value"] = 392245678915  # sdsf_seed
        DFN.params["boundaryFaces"]["value"] = [1, 1, 0, 0, 0, 0]
        DFN.params["keepOnlyLargestCluster"]["value"] = True
        # DFN.params['keepIsolatedFractures']['value'] = True
        # DFN.params['rkappa']['value']=0
        DFN.add_fracture_family(
            shape="rect",
            distribution="constant",
            probability=1,
            kappa=0.00001,
            theta=0,
            phi=0,
            constant=7.5 * (2**0.5) * 0.5,
            hy_variable="permeability",
            hy_function="correlated",
            hy_params={"alpha": 1e-13, "beta": 0.2},
        )

        start_time = time.time()
        DFN.make_working_directory(delete=True)
        DFN.check_input()
        DFN.print_domain_parameters()
        DFN.set_flow_solver("PFLOTRAN")
        DFN.create_network()
        DFN.output_report()
        elapsed_time = time.time() - start_time
        DFNGen_Time[Index_iiii, j] = elapsed_time

        start_time = time.time()
        DFN.mesh_network()
        elapsed_time = time.time() - start_time
        DFNMesh_Time[Index_iiii, j] = elapsed_time

        start_time = time.time()
        DFN.dfn_flow()
        elapsed_time = time.time() - start_time
        DFNFlow_Time[Index_iiii, j] = elapsed_time

        # DFN.dfn_trans()
    Index_iiii = Index_iiii + 1

    elapsed_time_iiii = time.time() - start_time_iiiii
    DFN_mesh_flow_total_time[Index_iiii-1, 0] = elapsed_time_iiii

file_name = src_path + "/ElapseTimeDFNWORKS.h5"

with h5py.File(file_name, "w") as file:
    # Create a dataset in the HDF5 file and write the NumPy array to it
    file.create_dataset("DFNGen_Time", data=DFNGen_Time)
    file.create_dataset("DFNMesh_Time", data=DFNMesh_Time)
    file.create_dataset("DFNFlow_Time", data=DFNFlow_Time)
    file.create_dataset("DFN_mesh_flow_total_time", data=DFN_mesh_flow_total_time)