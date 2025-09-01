import h5py
import numpy as np
from mayavi import mlab as ML
from numpy import linalg as LA
f = h5py.File('Flow_show.h5')
coordinate_3D = np.array(f['coordinate_3D'][:])
element_3D = np.array(f['element_3D'][:], dtype=int)
InletP = np.array(f['InletP'][:])
OutletP = np.array(f['OutletP'][:])
L_m = f['L_m'][:][0]
DomainDimensionRatio = f['DomainDimensionRatio'][:]
pressure_eles = np.array(f['pressure_eles'][:])
velocity_center_grid = np.array(f['velocity_center_grid'][:])
vec_norm = LA.norm(velocity_center_grid, axis=0)
f.close()
mesh = ML.triangular_mesh(coordinate_3D[0, :], coordinate_3D[1, :], coordinate_3D[2, :], np.transpose(element_3D-1), representation='wireframe', color=(0, 0, 0), line_width=1.0)
mesh.mlab_source.dataset.cell_data.scalars = np.transpose(pressure_eles)
mesh.mlab_source.dataset.cell_data.scalars.name = 'Cell data'
mesh.mlab_source.update()
mesh.parent.update()
mesh2 = ML.pipeline.set_active_attribute(mesh, cell_scalars='Cell data')
s2 = ML.pipeline.surface(mesh2, colormap='rainbow', opacity=0.8)
ML.outline(extent=[-0.5 * DomainDimensionRatio[0] * L_m, 0.5 * DomainDimensionRatio[0] * L_m, -0.5 * DomainDimensionRatio[1] * L_m, 0.5 * DomainDimensionRatio[1] * L_m, -0.5 * DomainDimensionRatio[2] * L_m, 0.5 * DomainDimensionRatio[2] * L_m])
ML.axes()
s2.module_manager.scalar_lut_manager.data_range = np.array([OutletP[0], InletP[0]])
ML.colorbar(object=s2, orientation='vertical', title='Hydraulic head [L]')
ML.xlabel('x (m)')
ML.ylabel('y (m)')
ML.zlabel('z (m)')
CenterELE = np.zeros([3, element_3D.shape[1]])
CenterELE[0, :] = 1.0 / 3.0 * (coordinate_3D[0, element_3D[0, :]-1] + coordinate_3D[0, element_3D[1, :]-1] + coordinate_3D[0, element_3D[2, :]-1])
CenterELE[1, :] = 1.0 / 3.0 * (coordinate_3D[1, element_3D[0, :]-1] + coordinate_3D[1, element_3D[1, :]-1] + coordinate_3D[1, element_3D[2, :]-1])
CenterELE[2, :] = 1.0 / 3.0 * (coordinate_3D[2, element_3D[0, :]-1] + coordinate_3D[2, element_3D[1, :]-1] + coordinate_3D[2, element_3D[2, :]-1])
ML.quiver3d(CenterELE[0, :], CenterELE[1, :], CenterELE[2, :], velocity_center_grid[0, :], velocity_center_grid[1, :], velocity_center_grid[2, :])
ML.show()
