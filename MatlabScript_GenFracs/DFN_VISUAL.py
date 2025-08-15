import h5py
import numpy as np
from mayavi import mlab as ML

f = h5py.File('DFN_VISUAL.h5')

verts = np.array(f['verts'][:])
Frac_NUM_verts = np.array(f['Frac_NUM_verts'][:])
verts = np.array(f['verts'][:])
L_m = f['L_m'][0]
DomainDimensionRatio = f['DomainDimensionRatio'][:]
NumClusters=np.array(f['NumClusters'][:])
ListClusters=[]
for i in range(NumClusters[0]):
	ListClusters.append(np.array(f["Cluster_" + str(i + 1)][:]))
PercolationClusters = np.array(f['PercolationClusters'][:])
Polar_Orientation = np.array(f['Polar_Orientation'][:])


intersections = np.array(f['intersections'][:])
NumIntersections = intersections.shape[1] if intersections.ndim != 1 else 1
if intersections.ndim == 1:
	intersections=intersections[:, np.newaxis]
Intersection =  np.concatenate((np.transpose(intersections[[0, 1, 2], :]), np.transpose(intersections[[3, 4, 5], :])), axis=0)
Connection_ = list()
Connection_ = [(Connection_ + [i, i + NumIntersections]) for i in range(NumIntersections)]
src = ML.pipeline.scalar_scatter(Intersection[:, 0], Intersection[:, 1], Intersection[:, 2])
src.mlab_source.dataset.lines = Connection_
src.update()
lines = ML.pipeline.stripper(src)
ML.pipeline.surface(lines, color=(1, 0, 0), line_width=5, opacity=1)

f.close()

NUM_Fracs = Frac_NUM_verts.shape[1]

poo = 0
structure_ = list()
scalar_t = np.zeros([verts.shape[1]])
for i in range(NUM_Fracs):
	NumVerticesSingleFrac = int(Frac_NUM_verts[0, i])
	pos_t = 0
	index_yt=0
	for j in range(NumClusters[0]):
		index_yt = np.where(ListClusters[j] == i + 1)
		if index_yt[0].size:
			pos_t = j
			break
	if not (pos_t + 1 in PercolationClusters):
		scalar_t[[s for s in range(poo, poo + NumVerticesSingleFrac)]] = pos_t
	else:
		scalar_t[[s for s in range(poo, poo + NumVerticesSingleFrac)]] = 1.5 * NumClusters
	for j in range(NumVerticesSingleFrac - 2):
		structure_ = structure_ + [[poo, poo + (j + 1) % NumVerticesSingleFrac, poo + (j + 2) % NumVerticesSingleFrac]]
	poo += NumVerticesSingleFrac

ML.triangular_mesh(verts[0, :], verts[1, :], verts[2, :], structure_, scalars=verts[2, :], opacity=1)
ML.outline(extent=[-0.5 * DomainDimensionRatio[0] * L_m, 0.5 * DomainDimensionRatio[0] * L_m, -0.5 * DomainDimensionRatio[1] * L_m, 0.5 * DomainDimensionRatio[1] * L_m, -0.5 * DomainDimensionRatio[2] * L_m, 0.5 * DomainDimensionRatio[2] * L_m])
ML.axes()
ML.colorbar(orientation='vertical')
ML.xlabel('x (m)')
ML.ylabel('y (m)')
ML.zlabel('z (m)')
ML.show()

ML.triangular_mesh(verts[0, :], verts[1, :], verts[2, :], structure_, scalars=scalar_t, opacity=1)
ML.outline(extent=[-0.5 * DomainDimensionRatio[0] * L_m, 0.5 * DomainDimensionRatio[0] * L_m, -0.5 * DomainDimensionRatio[1] * L_m, 0.5 * DomainDimensionRatio[1] * L_m, -0.5 * DomainDimensionRatio[2] * L_m, 0.5 * DomainDimensionRatio[2] * L_m])
ML.axes()
ML.colorbar(orientation='vertical')
ML.xlabel('x (m)')
ML.ylabel('y (m)')
ML.zlabel('z (m)')
ML.show()

import matplotlib.pyplot as plt
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(Polar_Orientation[0, :], Polar_Orientation[1, :], 'ko')
ax.set_xticks([0, 1.5708, 3.1416, 4.7124])
ax.set_xticklabels([r'0', r"$0.5\pi$", r"$\pi$", r"$1.5\pi$"])
ax.set_rmax(1.5708)
ax.set_rticks([0.5236, 1.0472, 1.5708])
ax.set_yticklabels([r"$0.25\pi$", r"$0.75\pi$", r"$\pi$"])
ax.set_rlabel_position(10)
ax.grid(True)
ax.set_title("Orientations", va='bottom')
plt.show()
