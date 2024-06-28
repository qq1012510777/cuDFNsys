from pydfnworks import *
import os
import numpy as np
import h5py
import time
# import random
src_path = os.getcwd()
dfnFlow_file = src_path + "/dfn_explicit.in"
dfnTrans_file = src_path + "/PTDFN_control.dat"
jobname = src_path + "/output"
DFN = DFNWORKS(
    jobname, dfnFlow_file=dfnFlow_file, dfnTrans_file=dfnTrans_file, ncpu=10
)
# sdsf_seed = int(time.time()) ^ int.from_bytes(os.urandom(4), byteorder='big')
DFN.params["domainSize"]["value"] = [30, 30, 30]
DFN.params["domainSizeIncrease"]["value"] = [0, 0, 0]
DFN.params["h"]["value"] = 0.1
DFN.params["stopCondition"]["value"] = 0
DFN.params["nPoly"]["value"] = 90
DFN.params["seed"]["value"] = 392245678915  # sdsf_seed
DFN.params["boundaryFaces"]["value"] = [0, 0, 1, 1, 0, 0]
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
DFN.dfn_gen()
DFN.mesh_network()
DFN.dfn_flow()
DFN.dfn_trans()



