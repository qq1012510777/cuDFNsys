import h5py
import os
import numpy as np
import sys

current_directory = os.getcwd()

if len(sys.argv) != 2:
    print("Please input the Percolation direction")
    sys.exit(1)
dir = int(sys.argv[1])
Axis = np.array(["X", "Y", "Z"])

q = np.zeros((3, 1))
HeadGradient = 0

f1 = h5py.File("./Class_DFN.h5")

Lm = np.array(f1["L"][0])

f1.close()

Num_Fracs = 0

P32_LargestCluster = 0
P32_Total = 0

P33_LargestCluster = 0
P33_Total = 0

Percolation_Status = 0
Permeability_Apparent = 0

HeadGradient = 1

f3 = h5py.File("./Class_DFN_original.h5")
Num_Fracs = np.array(f3["NumFractures"][0])

Area_and_aperture = np.zeros((Num_Fracs, 2))
for l in range(1, Num_Fracs + 1):
    Radius = np.array(f3["Fracture_" + str(l) + "/Radius"][0])
    Aperture = (np.array(f3["Fracture_" + str(l) + "/Conductivity"][0]) * 12) ** (
        1.0 / 3.0
    )
    Area_and_aperture[l - 1, 0] = (2**0.5 * Radius) ** 2
    Area_and_aperture[l - 1, 1] = Aperture
    # print(np.shape(Aperture))
    # print(np.shape(Area_and_aperture))
    # print("------------")

P32_Total = np.sum(Area_and_aperture[:, 0]) / (Lm**3)
P33_Total = np.sum(Area_and_aperture[:, 0] * Area_and_aperture[:, 1]) / (Lm**3)

NumPercolationCluster = f3["PercolationClusters"].size

if NumPercolationCluster != 0:
    Cluster_p = np.array([], dtype=int)
    for i in range(1, NumPercolationCluster + 1):
        PercolationClu = np.array(f3["PercolationClusters"][i - 1]) + 1
        Cluster_p_t = np.array(f3["Cluster_" + str(PercolationClu)][:])
        Cluster_p = np.concatenate((Cluster_p, Cluster_p_t), axis=0)
    P32_LargestCluster = np.sum(Area_and_aperture[Cluster_p, 0]) / (Lm**3)
    P33_LargestCluster = np.sum(
        Area_and_aperture[Cluster_p, 0] * Area_and_aperture[Cluster_p, 1]
    ) / (Lm**3)
else:
    P32_LargestCluster = 0
    P33_LargestCluster = 0
f3.close()


if not os.path.exists("./" + Axis[dir] + "_DarcyFlow_Finished"):
    print(current_directory, " is not finished~")
if not os.path.exists("./DFN_FLOW_VISUAL.h5"):
    # q[k - 1, k - 1] = 0
    Percolation_Status = 0
    Permeability_Apparent = 0
else:
    f2 = h5py.File("./DFN_FLOW_VISUAL.h5")
    velocity_center_grid = np.array(f2["velocity_center_grid"][:])
    ElementAperture = np.array(f2["ElementAperture"][:])
    f2.close()
    f3 = h5py.File("./Class_MESH.h5")
    coordinate_2D = np.array(f3["/group_mesh/Coordinate2D"][:])
    f3.close()
    f4 = h5py.File("./Class_FLOW.h5")
    Percolation_Status = 1
    Permeability_Apparent = np.array(f4["Permeability"][0])
    f4.close()
    L1 = np.linalg.norm(coordinate_2D[[1, 4], :] - coordinate_2D[[0, 3], :], axis=0)
    L2 = np.linalg.norm(coordinate_2D[[2, 5], :] - coordinate_2D[[1, 4], :], axis=0)
    L3 = np.linalg.norm(coordinate_2D[[0, 3], :] - coordinate_2D[[2, 5], :], axis=0)
    p = (L1 + L2 + L3) / 2
    Area = (p * (p - L1) * (p - L2) * (p - L3)) ** 0.5
    q_vector = velocity_center_grid * Area
    q[0, 0] = np.sum(q_vector[0, :]) / (Lm**3)
    q[1, 0] = np.sum(q_vector[1, :]) / (Lm**3)
    q[2, 0] = np.sum(q_vector[2, :]) / (Lm**3)
    # print(Num_Fracs)
    # print(P32_LargestCluster)
    # print(P32_Total)
    # print(P33_LargestCluster)
    # print(P33_Total)
    # print(Percolation_Status)
    # print(Permeability_Apparent)
    # print(i, j, k)
    # exit(1)
# if (np.sum(Percolation_Status) == 3):
#     print(q)
#     print(HeadGradient)
#     exit(1)
dataFIle = current_directory + "/data.h5"
lo = "w"
with h5py.File(dataFIle, lo) as hdf:
    hdf.create_dataset("Num_Fracs", data=Num_Fracs)
    hdf.create_dataset("P32_LargestCluster", data=P32_LargestCluster)
    hdf.create_dataset("P32_Total", data=P32_Total)
    hdf.create_dataset("P33_LargestCluster", data=P33_LargestCluster)
    hdf.create_dataset("P33_Total", data=P33_Total)
    hdf.create_dataset("Percolation_Status", data=Percolation_Status)
    hdf.create_dataset("Permeability_Apparent", data=Permeability_Apparent)
    hdf.create_dataset("q", data=q)
    hdf.create_dataset("Lm", data=Lm)
    hdf.create_dataset("HeadGradient", data=HeadGradient)
# if j == 20:
#    print(Permeability_j)
#    print(q)
#    print(np.linalg.inv(HeadGradient))
#    exit(1)
