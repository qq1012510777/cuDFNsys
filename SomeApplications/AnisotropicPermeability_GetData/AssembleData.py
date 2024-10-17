import h5py
import os
import numpy as np
import sys

# Check if enough arguments are passed
if len(sys.argv) != 3:
    print("Usage: python script_name.py <start> <end>")
    print("Please input the start and end set number")
    sys.exit(1)

# Get the start and end values from command-line arguments
start = int(sys.argv[1])
end = int(sys.argv[2])

# Get the current working directory
current_directory = os.getcwd()

try:
    os.mkdir(current_directory + "/data")
    print(f"Directory '{current_directory + "/data"}' created.")
except FileExistsError:
    print(f"Directory '{current_directory + "/data"}' already exists.")

dataFIle_assem = current_directory + "/data/data.h5"
SetNO = np.arange(start, end + 1)

for i in SetNO:
    print("addressing set_" + str(i))
    for j in np.arange(1, 31):
        Axis = np.array(["X", "Y", "Z"])

        for k in range(1, 4):
            FilenameDFN = (
                current_directory
                + "/set"
                + str(i).zfill(5)
                + "/DFN_"
                + Axis[k - 1]
                + "_"
                + str(j).zfill(3)
            )

            try:
                dataFIle = FilenameDFN + "/data.h5"

                f1 = h5py.File(dataFIle)
                HeadGradient = np.array(f1['HeadGradient'])
                Lm = np.array(f1['Lm'])
                Num_Fracs = np.array(f1['Num_Fracs'])
                P32_LargestCluster = np.array(f1['P32_LargestCluster'])
                P32_Total  = np.array(f1['P32_Total'])
                P33_LargestCluster = np.array(f1['P33_LargestCluster'])
                P33_Total = np.array(f1['P33_Total'])
                Percolation_Status = np.array(f1['Percolation_Status'])
                Permeability_Apparent = np.array(f1['Permeability_Apparent'])
                q = np.array(f1['q'][:])
                f1.close()


                lo = 'w' if i == start and j == 1 and k == 1 else 'a'

                with h5py.File(dataFIle_assem, lo) as hdf:
                    group1 = hdf.create_group("/set"
                        + str(i).zfill(5)
                        + "_DFN_" + Axis[k - 1] + "_"
                        + str(j).zfill(3))
                    group1.create_dataset('HeadGradient', data=HeadGradient)
                    group1.create_dataset('Lm', data=Lm)
                    group1.create_dataset('Num_Fracs', data=Num_Fracs)
                    group1.create_dataset('P32_LargestCluster', data=P32_LargestCluster)
                    group1.create_dataset('P32_Total', data=P32_Total)
                    group1.create_dataset('P33_LargestCluster', data=P33_LargestCluster)
                    group1.create_dataset('P33_Total', data=P33_Total)
                    group1.create_dataset('Percolation_Status', data=Percolation_Status)
                    group1.create_dataset('Permeability_Apparent', data=Permeability_Apparent)
                    group1.create_dataset('q', data=q)
            except Exception as e:
                print("!!! ~~~", FilenameDFN, " is wrong ~")
                print(f"An error occurred: {e}")
                