clc;
clear all
close all

load('MHFEM_1.mat');
L = 0.5 * 100;
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

5082
5081
5085
4974
4922
5005
4834
4850
4970
5008
5051
4999
    ];
patch('Vertices', coordinate_3D, 'Faces', element_3D(Tri, :), 'FaceVertexCData', pressure_eles(Tri, :), 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 0.9); view(3); hold on
colorbar; hold on
xlim([5.60, 5.68]);
ylim([-0.633, -0.58]);
zlim([-37.79, -37.69]);
hold on

Trajectory = [	5.6640739325474456222764274571090936660767, -0.5890290414717043177006416954100131988525, -37.6914262245349434010677214246243238449097
	5.6027291215561731974048598203808069229126, -0.6309092559992564019921701401472091674805, -37.7897802193474134924144891556352376937866

    ];
plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-r', 'linewidth', 2); hold on

% scatter3(33.2435605634269464303542918059974908828735, -26.6327849633078450608536513755097985267639, 41.4994501729028613112859602551907300949097, ...
%     '^', 'k');

Trajectory = [		5.6640545454506110445436206646263599395752, -0.5890422770781214012458804063498973846436, -37.6914573078228372082776331808418035507202
	5.6027291215561731974048598203808069229126, -0.6309092559992564019921701401472091674805, -37.7897802193474134924144891556352376937866
];
plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-b', 'linewidth', 2); hold on

% Trajectory = [	4.2162086765735278248712347703985869884491, 22.6895633559991054539750621188431978225708, -19.0073204054348074976132920710369944572449
% 	4.2198070351265766220194564084522426128387, 22.6767084168575401292855531210079789161682, -19.0147733131078879864617192652076482772827];
% plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-g', 'linewidth', 2); hold on
% 
% Trajectory = [	4.2177515988141305314229612122289836406708, 22.6840513506748138183866103645414113998413, -19.0105161005931506679189624264836311340332
% 	4.2198070351265766220194564084522426128387, 22.6767084168575401292855531210079789161682, -19.0147733131078879864617192652076482772827];
% plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-c', 'linewidth', 2); hold on
% 
% Trajectory = [		4.2186596170619878876095754094421863555908, 22.6808075051425888091216620523482561111450, -19.0123967848181898432358138961717486381531
% 	4.2198070351265766220194564084522426128387, 22.6767084168575401292855531210079789161682, -19.0147733131078879864617192652076482772827];
% plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-k', 'linewidth', 2); hold on