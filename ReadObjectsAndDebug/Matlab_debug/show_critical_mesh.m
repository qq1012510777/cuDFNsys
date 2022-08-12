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
1697
1602
2975

    ];
patch('Vertices', coordinate_3D, 'Faces', element_3D(Tri, :), 'FaceVertexCData', pressure_eles(Tri, :), 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 0.9); view(3); hold on
colorbar; hold on
% xlim([-10.23, -10.21]);
% ylim([-0.58, -0.56]);
% zlim([9.64, 9.66]);
hold on

Trajectory = [	-10.2230030211751810043097066227346658706665, -0.5771106071272665838378657099383417516947, 9.6509090655413363180059604928828775882721
	-10.2255435687823208468216762412339448928833, -0.5782937448162919213956456587766297161579, 9.6465624918257919517827758681960403919220
    ];
plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-r', 'linewidth', 2); hold on

Trajectory = [			-10.2234376270102096384562173625454306602478, -0.5773130038632654503771846066229045391083, 9.6501655068322325092822211445309221744537
	-10.2204196476527098269571069977246224880219, -0.5743292065775564836371813726145774126053, 9.6495594416076606592014286434277892112732
    ];
plot3(Trajectory(:, 1), Trajectory(:, 2), Trajectory(:, 3), '*-b', 'linewidth', 2); hold on

