import sys
sys.path.append("/home/torresLab/yintingchang/conda/ENTER")
import h5py
import numpy as np

f1=h5py.File('ParticlePositionResult/DispersionInfo.h5', 'r+')
NumSteps=np.array(f1['NumOfSteps'][0])
f1.close()

f2=h5py.File('Dispersion_MeanSquareDisplacement.h5', 'r+')
io=1
try:
    while io > 0:
        NumSteps=NumSteps+1
        # print(('Mean_MSD_cMSD_'+str(NumSteps).zfill(10)))
        del f2[('Mean_MSD_cMSD_'+str(NumSteps).zfill(10))]

except:
    print("Finished")

f2.close()