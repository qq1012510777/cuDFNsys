clc;
close all;
clear all;
currentPath = fileparts(mfilename('fullpath'));
coordinate_3D = h5read([currentPath, '/Mesh_show.h5'], '/coordinate_3D');
element_3D = h5read([currentPath, '/Mesh_show.h5'], '/element_3D');
L = h5read([currentPath, '/Mesh_show.h5'], '/L_m');
DomainDimensionRatio = h5read([currentPath, '/Mesh_show.h5'], '/DomainDimensionRatio');
cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, L, -L; -L,L, L];
cube_frame(:, 1) = 0.5 .* cube_frame(:, 1) .* DomainDimensionRatio(1); cube_frame(:, 2) = 0.5 .* cube_frame(:, 2) .* DomainDimensionRatio(2); cube_frame(:, 3) = 0.5 .* cube_frame(:, 3) .* DomainDimensionRatio(3);
figure(1); view(3); title('DFN mesh'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on
patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on
patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', coordinate_3D(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1); hold on
axis([-1.1 / 2 * DomainDimensionRatio(1) * L,  1.1 / 2 * DomainDimensionRatio(1) * L, -1.1 / 2 * DomainDimensionRatio(2) * L, 1.1 / 2 * DomainDimensionRatio(2) * L, -1.1 / 2 * DomainDimensionRatio(3) * L, 1.1 / 2 * DomainDimensionRatio(3) * L]);
pbaspect([DomainDimensionRatio]); hold on
