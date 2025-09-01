import h5py
import numpy as np
from mayavi import mlab as ML
f = h5py.File('Mesh_show.h5')
coordinate_3D = np.array(f['coordinate_3D'][:])
element_3D = np.array(f['element_3D'][:])
L_m = f['L_m'][:][0]
DomainDimensionRatio = f['DomainDimensionRatio'][:]
f.close()
ML.triangular_mesh(coordinate_3D[0, :], coordinate_3D[1, :], coordinate_3D[2, :], np.transpose(element_3D - 1), scalars=coordinate_3D[2, :], opacity=0.8)
ML.triangular_mesh(coordinate_3D[0, :], coordinate_3D[1, :], coordinate_3D[2, :], np.transpose(element_3D-1), representation='wireframe', color=(0, 0, 0), line_width=1.0)
ML.outline(extent=[-0.5 * DomainDimensionRatio[0] * L_m, 0.5 * DomainDimensionRatio[0] * L_m, -0.5 * DomainDimensionRatio[1] * L_m, 0.5 * DomainDimensionRatio[1] * L_m, -0.5 * DomainDimensionRatio[2] * L_m, 0.5 * DomainDimensionRatio[2] * L_m])
ML.axes()
ML.colorbar(orientation='vertical')
ML.xlabel('x (m)')
ML.ylabel('y (m)')
ML.zlabel('z (m)')
ML.show()
