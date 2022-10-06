import h5py
import mpl_toolkits.mplot3d as a3
import matplotlib
import matplotlib.colors as colors
import pylab as pl
import numpy as np
from itertools import product, combinations

f = h5py.File('DFN_II.h5')
for key in f.keys():
    print(f[key].name)

verts = np.array(f['verts'][:])
Frac_NUM_verts = np.array(f['Frac_NUM_verts'][:])
ListClusters = np.array(f['ListClusters'][:])
PercolationClusters = np.array(f['PercolationClusters'][:])
Polar_Orientation = np.array(f['Polar_Orientation'][:])
Rs_ = np.array(f['R'][:])
intersections = np.array(f['intersections'][:])
verts = np.array(f['verts'][:])
f.close()

NUM_Fracs = Frac_NUM_verts.shape[1]

init_t = 0

fig = pl.figure()
ax = a3.Axes3D(fig, auto_add_to_figure=False)
L = 30

minX = min(verts[0, :])
minY = min(verts[1, :])
minZ = min(verts[2, :])
maxX = max(verts[0, :])
maxY = max(verts[1, :])
maxZ = max(verts[2, :])
minO = min([minX, minY, minZ]) - 10
maxO = max([maxX, maxY, maxZ]) + 10

ax.axes.set_xlim3d(left=minO, right=maxO)
ax.axes.set_ylim3d(bottom=minO, top=maxO)
ax.axes.set_zlim3d(bottom=minO, top=maxO)
ax.axes.set_xlabel('$x$ (m)', fontweight='bold')
ax.axes.set_ylabel('$y$ (m)', fontweight='bold')
ax.axes.set_zlabel('$z$ (m)', fontweight='bold')
fig.suptitle('Discrete fracture network\n',
             fontsize = 14, fontweight ='bold')

# intersection
fig.add_axes(ax)
for i in range(intersections.shape[1]):
    vtx1 = np.zeros([2, 3]);
    vtx1[0, :] = np.transpose(np.array(intersections[[0, 1, 2], i]))
    vtx1[1, :] = np.transpose(np.array(intersections[[3, 4, 5], i]))
    Intersection_ = a3.art3d.Line3DCollection([vtx1], linewidths=5)
    Intersection_.set_color(colors.rgb2hex([0, 0, 1]))
    Intersection_.set_edgecolor('r')
    ax.add_collection3d(Intersection_)

fig.add_axes(ax)
lambdas = range(1, 9)
for i in range(NUM_Fracs):
    vtx = np.array(verts[:, init_t:init_t + int(Frac_NUM_verts[0, i])])
    # print(np.transpose(vtx))
    init_t += int(Frac_NUM_verts[0, i])
    Frac = a3.art3d.Poly3DCollection([np.transpose(vtx)], alpha=1)
    Frac.set_color(colors.rgb2hex(np.random.rand(3)))
    Frac.set_edgecolor('k')
    ax.add_collection3d(Frac, zs=lambdas, zdir='y')


# the domain
# r = [-0.5 * L, 0.5 * L]
# for s, e in combinations(np.array(list(product(r, r, r))), 2):
#     if np.sum(np.abs(s - e)) == r[1] - r[0]:
#         ax.plot3D(*zip(s, e), color="black")

pl.show()
