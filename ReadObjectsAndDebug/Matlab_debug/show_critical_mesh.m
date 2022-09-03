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
48265
52002
    ];
patch('Vertices', coordinate_3D, 'Faces', element_3D(Tri, :), 'FaceVertexCData', pressure_eles(Tri, :), 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 0.9); view(3); hold on
colorbar; hold on
xlim([37.800 37.809]);
ylim([-15.1568 -15.15608]);
zlim([14.280 14.288]);
hold on

Trajectory = [	37.8002693876271393946808530017733573913574, -15.1568026486496769678069540532305836677551, 14.2876537276872941362171332002617418766022
	37.8088966932151464561684406362473964691162, -15.1560986184463946102596310083754360675812, 14.2802897586648143146703660022467374801636
    ];
plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-r', 'linewidth', 2); hold on

% plot3(coordinate_3D([7819, 927], 1), coordinate_3D([7819, 927], 2), coordinate_3D([7819, 927], 3), '*-k', 'linewidth', 2); hold on


% scatter3(33.2435605634269464303542918059974908828735, -26.6327849633078450608536513755097985267639, 41.4994501729028613112859602551907300949097, ...
%     '^', 'k');

% Trajectory = [24.4275466154353608771998551674187183380127, 2.4310887384975732317116126068867743015289, 21.2483743618818543552606570301577448844910
% 	24.4442131104148039355550281470641493797302, 2.4520563272697089018947735894471406936646, 21.2143054842207021692956914193928241729736
% ];
% plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), 'o-b', 'linewidth', 2); hold on
% 
% Trajectory = [	24.4279640857514550589257851243019104003906, 2.4316139446625948394853367062751203775406, 21.2475209883831475110582687193527817726135
% 	24.3892219788303883376556768780574202537537, 2.4161427977722977900043588306289166212082, 21.2408214144285238944576121866703033447266];
% 
% plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-g', 'linewidth', 2); hold on
% % 
% Trajectory = [	4.2177515988141305314229612122289836406708, 22.6840513506748138183866103645414113998413, -19.0105161005931506679189624264836311340332
% 	4.2198070351265766220194564084522426128387, 22.6767084168575401292855531210079789161682, -19.0147733131078879864617192652076482772827];
% plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-c', 'linewidth', 2); hold on
% 
% Trajectory = [		4.2186596170619878876095754094421863555908, 22.6808075051425888091216620523482561111450, -19.0123967848181898432358138961717486381531
% 	4.2198070351265766220194564084522426128387, 22.6767084168575401292855531210079789161682, -19.0147733131078879864617192652076482772827];
% plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-k', 'linewidth', 2); hold on