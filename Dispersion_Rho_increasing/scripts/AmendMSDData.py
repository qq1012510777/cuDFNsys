import sys
sys.path.append("/home/torresLab/yintingchang/conda/ENTER")
import h5py
import numpy as np

f1=h5py.File('ParticlePositionResult/DispersionInfo.h5', 'r+')
NumSteps=np.array(f1['NumOfSteps'][:])
f1.close()

f2=h5py.File('Dispersion_MeanSquareDisplacement.h5', 'r+')
f3=h5py.File('new_Dispersion_MeanSquareDisplacement.h5', "w")    # mode = {'w', 'r', 'a'}
for i in range(0,NumSteps[0]+1):
    UI=np.array(f2[('Mean_MSD_cMSD_'+str(i).zfill(10))][:])
    d = f3.create_dataset(('Mean_MSD_cMSD_'+str(i).zfill(10)), data=UI)

f2.close()
f3.close()
