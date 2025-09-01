clc;
close all;
clear all;
currentPath = fileparts(mfilename('fullpath'));
coordinate_3D = h5read([currentPath, '/Flow_show.h5'], '/coordinate_3D');
element_3D = h5read([currentPath, '/Flow_show.h5'], '/element_3D');
velocity_center_grid = h5read([currentPath, '/Flow_show.h5'], '/velocity_center_grid');
pressure_eles = h5read([currentPath, '/Flow_show.h5'], '/pressure_eles');
VelocityNorm = vecnorm(velocity_center_grid');

L = h5read([currentPath, '/Flow_show.h5'], '/L_m');
DomainDimensionRatio = h5read([currentPath, '/Flow_show.h5'], '/DomainDimensionRatio');
cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, L, -L; -L,L, L];
cube_frame(:, 1) = 0.5 .* cube_frame(:, 1) .* DomainDimensionRatio(1); cube_frame(:, 2) = 0.5 .* cube_frame(:, 2) .* DomainDimensionRatio(2); cube_frame(:, 3) = 0.5 .* cube_frame(:, 3) .* DomainDimensionRatio(3);
figure(1); view(3); title('DFN flow (mhfem)'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on
patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on

axis([-1.1 / 2 * DomainDimensionRatio(1) * L,  1.1 / 2 * DomainDimensionRatio(1) * L, -1.1 / 2 * DomainDimensionRatio(2) * L, 1.1 / 2 * DomainDimensionRatio(2) * L, -1.1 / 2 * DomainDimensionRatio(3) * L, 1.1 / 2 * DomainDimensionRatio(3) * L]);
pbaspect([DomainDimensionRatio]); hold on
patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', pressure_eles, 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 1); view(3); hold on
Cb = colorbar;
Cb.Title.String = 'Hydraulic head [L]';
caxis([0, 1]);

CenterELE = zeros(size(element_3D, 1), 3);
CenterELE(:, 1) = 1/3 * (coordinate_3D(element_3D(:, 1), 1) + coordinate_3D(element_3D(:, 2), 1) + coordinate_3D(element_3D(:, 3), 1));
CenterELE(:, 2) = 1/3 * (coordinate_3D(element_3D(:, 1), 2) + coordinate_3D(element_3D(:, 2), 2) + coordinate_3D(element_3D(:, 3), 2));
CenterELE(:, 3) = 1/3 * (coordinate_3D(element_3D(:, 1), 3) + coordinate_3D(element_3D(:, 2), 3) + coordinate_3D(element_3D(:, 3), 3));
quiver3(CenterELE(:, 1), CenterELE(:, 2), CenterELE(:, 3), velocity_center_grid(:, 1),velocity_center_grid(:, 2),velocity_center_grid(:, 3), 4, 'LineWidth', 1.5, 'color', 'r');
ElementAperture = h5read([currentPath, '/Flow_show.h5'], '/ElementAperture');
VelocityNorm = [vecnorm(velocity_center_grid')]' ./ ElementAperture;
len1=[vecnorm([coordinate_3D(element_3D(:, 1), :) - coordinate_3D(element_3D(:, 2), :)]')]';
len2=[vecnorm([coordinate_3D(element_3D(:, 3), :) - coordinate_3D(element_3D(:, 2), :)]')]';
len3=[vecnorm([coordinate_3D(element_3D(:, 1), :) - coordinate_3D(element_3D(:, 3), :)]')]';
P_ss = (len1+len2+len3)*0.5;
Area_ss=(P_ss .* (P_ss-len1) .* (P_ss-len2) .* (P_ss-len3)) .^ 0.5;
meanFractureVelocity = sum(VelocityNorm .* Area_ss .* ElementAperture) ./ (sum(Area_ss .* ElementAperture))
