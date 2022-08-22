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

18363,
13753,
13581,
13655
    ];
patch('Vertices', coordinate_3D, 'Faces', element_3D(Tri, :), 'FaceVertexCData', pressure_eles(Tri, :), 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 0.9); view(3); hold on
colorbar; hold on
xlim([32.16, 32.17]);
ylim([-37.13, -37.126]);
zlim([12.924, 12.93]);
hold on

Trajectory = [			32.1689254637597557007211435120552778244019, -37.1292192330497954344536992721259593963623, 12.9270378962602414674165629548951983451843
	32.1691057437138354657690797466784715652466, -37.1291383062955233640423102770000696182251, 12.9269140011489547248402232071384787559509

    ];
plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-r', 'linewidth', 2); hold on

% scatter3(33.2435605634269464303542918059974908828735, -26.6327849633078450608536513755097985267639, 41.4994501729028613112859602551907300949097, ...
%     '^', 'k');

Trajectory = [	32.1690673504379844871436944231390953063965, -37.1291555408417366379580926150083541870117, 12.9269403864406839943512750323861837387085
	32.1690650809638469809215166606009006500244, -37.1292013892276528963520831894129514694214, 12.9269214119359201475845111417584121227264
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