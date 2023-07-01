clc;
close all;
clear all;
currentPath = fileparts(mfilename('fullpath'));
coordinate_3D = h5read([currentPath, '/MHFEM_1.h5'], '/coordinate_3D');
element_3D = h5read([currentPath, '/MHFEM_1.h5'], '/element_3D');
velocity_center_grid = h5read([currentPath, '/MHFEM_1.h5'], '/velocity_center_grid');
pressure_eles = h5read([currentPath, '/MHFEM_1.h5'], '/pressure_eles');
L = h5read([currentPath, '/MHFEM_1.h5'], '/L_m');
DomainDimensionRatio = h5read([currentPath, '/MHFEM_1.h5'], '/DomainDimensionRatio');
cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, L, -L; -L,L, L];
cube_frame(:, 1) = 0.5 .* cube_frame(:, 1) .* DomainDimensionRatio(1); cube_frame(:, 2) = 0.5 .* cube_frame(:, 2) .* DomainDimensionRatio(2); cube_frame(:, 3) = 0.5 .* cube_frame(:, 3) .* DomainDimensionRatio(3);
figure(1); view(3); title('DFN flow (mhfem)'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on
patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on

axis([-1.1 / 2 * DomainDimensionRatio(1) * L,  1.1 / 2 * DomainDimensionRatio(1) * L, -1.1 / 2 * DomainDimensionRatio(2) * L, 1.1 / 2 * DomainDimensionRatio(2) * L, -1.1 / 2 * DomainDimensionRatio(3) * L, 1.1 / 2 * DomainDimensionRatio(3) * L]);
pbaspect([DomainDimensionRatio]); hold on
patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', pressure_eles, 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 1); colorbar; view(3); hold on
caxis([79.999, 100]);

CenterELE = zeros(size(element_3D, 1), 3);
CenterELE(:, 1) = 1/3 * (coordinate_3D(element_3D(:, 1), 1) + coordinate_3D(element_3D(:, 2), 1) + coordinate_3D(element_3D(:, 3), 1));
CenterELE(:, 2) = 1/3 * (coordinate_3D(element_3D(:, 1), 2) + coordinate_3D(element_3D(:, 2), 2) + coordinate_3D(element_3D(:, 3), 2));
CenterELE(:, 3) = 1/3 * (coordinate_3D(element_3D(:, 1), 3) + coordinate_3D(element_3D(:, 2), 3) + coordinate_3D(element_3D(:, 3), 3));
quiver3(CenterELE(:, 1), CenterELE(:, 2), CenterELE(:, 3), velocity_center_grid(:, 1),velocity_center_grid(:, 2),velocity_center_grid(:, 3), 4, 'LineWidth', 1.5, 'color', 'r');
VelocityNorm = vecnorm(velocity_center_grid');
ElementAperture = h5read([currentPath, '/MHFEM_1.h5'], '/ElementAperture');
VelocityNorm = [vecnorm(velocity_center_grid')]' ./ ElementAperture;
meanFractureVelocity = mean(VelocityNorm)
edge_attri = h5read([currentPath, '/MHFEM_1.h5'], '/edge_attri');
inlet_loc = find(edge_attri(:, 4)==0);
inlet_loc = [inlet_loc; find(edge_attri(:, 5)==0)];
inlet_loc = [inlet_loc; find(edge_attri(:, 6)==0)];
meanFractureVelocity_Inlet = mean(VelocityNorm(inlet_loc))

figure(2)
patch('Vertices', coordinate_3D, 'Faces', element_3D([1064, 224], :),......
    'FaceVertexCData', coordinate_3D(:, 3), 'FaceColor', ...
    'flat', 'EdgeAlpha', 0.9, 'facealpha', 0);
view(3); hold on

trajectory3D = [-11.3334378817644356729488208657130599021912, -0.0502831231025270136703042567205557134002, -0.1381517441859155259642477631132351234555
	-11.2767477207344377632125542731955647468567, 0.0001922872823749333830774949083419755880, 0.0005283049620983399632598298545360648859];
plot3(trajectory3D(:, 1), trajectory3D(:, 2), trajectory3D(:, 3), '-'); hold on

trajectory3D = [-11.2769636832579553953337381244637072086334, -0.0000000000000022083757922067703156254894, 0.0000000000000060674625736599419810237497
	-11.2771796457814676983844037749804556369781, -0.0001922872823801448176599809825049192114, 0.0005283049621126582624103251717428975098];
plot3(trajectory3D(:, 1), trajectory3D(:, 2), trajectory3D(:, 3), '-+'); hold on