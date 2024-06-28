from pydfnworks import *
import os
import numpy as np
import h5py
import time

import random

src_path = os.getcwd()

dfnFlow_file = src_path + "/dfn_explicit.in"
dfnTrans_file = src_path + "/PTDFN_control.dat"

NumFractures_init = 1500
FracIncreament = 200

Num_MCTimes = 30
Num_FracIncrements = 10
NumCPUThreads = 10

DFNGen_Time = np.zeros((1, Num_FracIncrements, Num_MCTimes))

Index_iiii = 0
for i in range(NumCPUThreads, NumCPUThreads + 1):
    for j in range(Num_FracIncrements):
        for k in range(Num_MCTimes):
            jobname = (
                src_path
                + "/output_CPUThreads_"
                + str(i + 1)
                + "_NumFracs_"
                + str(NumFractures_init + j * FracIncreament)
                + "_MC_NO_"
                + str(k + 1)
            )
            DFN = DFNWORKS(
                jobname,
                dfnFlow_file=dfnFlow_file,
                dfnTrans_file=dfnTrans_file,
                ncpu=i + 1,
            )
            sdsf_seed = (
                int(time.time()) ^ int.from_bytes(os.urandom(4), byteorder="big") + k
            )
            DFN.params["domainSize"]["value"] = [100, 100, 100]
            DFN.params["domainSizeIncrease"]["value"] = [0, 0, 0]
            DFN.params["h"]["value"] = 0.1
            DFN.params["stopCondition"]["value"] = 0
            DFN.params["nPoly"]["value"] = NumFractures_init + j * FracIncreament
            DFN.params["seed"]["value"] = sdsf_seed
            DFN.params["boundaryFaces"]["value"] = [1, 1, 0, 0, 0, 0]
            # DFN.params["keepOnlyLargestCluster"]["value"] = True
            DFN.params["keepIsolatedFractures"]["value"] = True
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

            DFN.make_working_directory(delete=True)
            DFN.check_input()
            DFN.print_domain_parameters()
            DFN.set_flow_solver("PFLOTRAN")
            start_time = time.time()
            DFN.create_network()
            elapsed_time = time.time() - start_time
            DFNGen_Time[Index_iiii, j, k] = elapsed_time

            # DFN.output_report()

            # DFN.dfn_trans()
    Index_iiii = Index_iiii + 1

file_name = src_path + "/ElapseTime_DFNGen_DFNWORKS.h5"

with h5py.File(file_name, "w") as file:
    # Create a dataset in the HDF5 file and write the NumPy array to it
    Index_iiii = 0
    for i in range(NumCPUThreads, NumCPUThreads+1):
        file.create_dataset(
            "DFNGen_Time_NumThreads_" + str(i + 1), data=DFNGen_Time[Index_iiii, :, :]
        )
        Index_iiii = Index_iiii + 1
