clc;
clear all
close all

load('MHFEM_1.mat');
L = 0.5 * 30;
cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; -L, L, -L; -L, L, L];
% figure(1); view(3); title('DFN flow (mhfem)'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on
% patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4; 5 6 7 8; 9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on
%
% xlim([-1.1 * L, 1.1 * L]);
% ylim([-1.1 * L, 1.1 * L]);
% zlim([-1.1 * L, 1.1 * L]);
% hold on
% patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', pressure_eles, 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 1); colorbar; view(3); hold on
% caxis([20, 100]);

CenterELE = zeros(size(element_3D, 1), 3);
CenterELE(:, 1) = 1/3 * (coordinate_3D(element_3D(:, 1), 1) + coordinate_3D(element_3D(:, 2), 1) + coordinate_3D(element_3D(:, 3), 1));
CenterELE(:, 2) = 1/3 * (coordinate_3D(element_3D(:, 1), 2) + coordinate_3D(element_3D(:, 2), 2) + coordinate_3D(element_3D(:, 3), 2));
CenterELE(:, 3) = 1/3 * (coordinate_3D(element_3D(:, 1), 3) + coordinate_3D(element_3D(:, 2), 3) + coordinate_3D(element_3D(:, 3), 3));
% quiver3(CenterELE(:, 1), CenterELE(:, 2), CenterELE(:, 3), velocity_center_grid(:, 1),velocity_center_grid(:, 2),velocity_center_grid(:, 3), 4, 'LineWidth', 1.5, 'color', 'r');

figure(2);
view(3)
Tri = [
    6982
    6993
    6984
    6082
    6315
    6956
    ];
patch('Vertices', coordinate_3D, 'Faces', element_3D(Tri, :), 'FaceVertexCData', pressure_eles(Tri, :), 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 0.9); view(3); hold on
colorbar; hold on
xlim([-10.1, -10.07]);
ylim([-3.58, -3.56]);
zlim([13.64, 13.67]);
hold on

Trajectory = [-10.0990895511251288496623601531609892845154, -3.5748917887509983337679386750096455216408, 13.6620018155939195736436886363662779331207
	-10.0981126192544348185720082256011664867401, -3.5743411631315091803173800144577398896217, 13.6578603526220714314831639057956635951996
];
plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-r', 'linewidth', 2); hold on

Trajectory = [		-10.0989790758293160877201444236561655998230, -3.5748295218416914487136182287940755486488, 13.6615334826828149772381948423571884632111
	-10.0981126192544348185720082256011664867401, -3.5743411631315091803173800144577398896217, 13.657860352622071431483163905795663595199];
plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-b', 'linewidth', 2); hold on

Trajectory = [		-10.0989142972907242068458799622021615505219, -3.5747930108791390324540770961903035640717, 13.6612588699621788634885888313874602317810
	-10.0981126192544348185720082256011664867401, -3.5743411631315091803173800144577398896217, 13.6578603526220714314831639057956635951996
];
plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-k', 'linewidth', 2); hold on

