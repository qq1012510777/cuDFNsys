import h5py
import numpy as np

file_name = 'DispersionInfo.h5'

f1 = h5py.File(file_name, 'r+')

Delta_T = f1['Delta_T']
Delta_T[...] = 200

BlockNOPresent = f1['BlockNOPresent']
BlockNOPresent[...] = 0

NumOfSteps = f1['NumOfSteps']
NumOfSteps[...] = 0

f1.close()