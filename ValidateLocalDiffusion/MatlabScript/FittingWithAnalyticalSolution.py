import matplotlib.pyplot as plt
import numpy as np
import pylab as mpl
from scipy import special
from scipy.optimize import curve_fit
import math
import h5py


def func(t, V, D, Lm):
    return 0.5 * (special.erfc(
        (Lm - V * t) /
        (2 * np.sqrt(D * t))) + np.exp(V * Lm / D) * special.erfc(
            (Lm + V * t) / (2 * np.sqrt(D * t))))


f1 = h5py.File('./ParticlePositionResult/ParticlePosition_WhichStepDoesTheParticleReached.h5')
f2 = h5py.File('./ParticlePositionResult/DispersionInfo.h5')

FPT = np.array(f1['WhichStepDoesTheParticleReached'][0, :])
minusOneIndices = np.where(FPT == -1)
FPT = np.delete(FPT, minusOneIndices)
FirstArrivalTime = int(np.min(FPT))

Dm = np.array(f2['Dispersion_local'])

DeltaT = np.array(f2['Delta_T'][0])

FPT = FPT * DeltaT
SizeFPT = np.shape(FPT)

counts, bin_edges = np.histogram(FPT, bins=50)

KS = (len(str(int(DeltaT)))) + 1

bin_edges[1:] = bin_edges[1:]
cdf = np.cumsum(counts) / (SizeFPT[0])

popt = 0
pcov = 0
found = 0
LoopTime=11
RSqTol = 0.88
for i in range(0, LoopTime):
    popt, pcov = curve_fit(lambda t, V, D: func(bin_edges[1:], V, D, 1),
                        bin_edges[1:], cdf,
                        [10**(-(KS + i)), 10**(-(KS + i))])
    residuals = cdf - func(bin_edges[1:], popt[0], popt[1], 1)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((cdf - np.mean(cdf))**2)
    r_squared = 1 - (ss_res / ss_tot)
    if (r_squared >= RSqTol):
        found = 1
        break               
if (found == 0):
    for i in range(-LoopTime, 0):
        # print(i, KS)
        popt, pcov = curve_fit(lambda t, V, D: func(bin_edges[1:], V, D, 1),
                            bin_edges[1:], cdf,
                            [10**(-(KS + i)), 10**(-(KS + i))])
        residuals = cdf - func(bin_edges[1:], popt[0], popt[1], 1)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((cdf - np.mean(cdf))**2)
        r_squared = 1 - (ss_res / ss_tot)
        if (r_squared >= RSqTol):
            break  
if (found == 0):
    exit(': : : Cannot find optimizing fitting')

plt.plot(bin_edges[1:], cdf, marker="o", linestyle='None')
plt.plot(bin_edges[1:], func(bin_edges[1:], popt[0], popt[1], 1), marker="v")
print('Dm: ' +  str(Dm) +  ', meanV = '  +  ', Vt = ' + str(popt[0]) + ', Dt = ' + str(popt[1]) + ', RSquare = ' + str(r_squared))

plt.xlabel('Time')
plt.ylabel('Frequency')

plt.show()
f1.close()
f2.close()