clc;
close all;
clear all;
load('particle.mat');
L = 0.5 * 30;
cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, L, -L; -L,L, L];
figure(1); view(3); title('DFN flow (mhfem)'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on
patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on

patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', pressure_eles, 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 1); colorbar; view(3); hold on
caxis([20, 100]);
xlim([-(0.1 * L + L), (0.1 * L + L)])
ylim([ -(0.1 * L + L), (0.1 * L + L) ])
zlim([ -(0.1 * L + L), (0.1 * L + L) ]);hold on

S = load('ParticlePosition.mat');
Matrx3D_pso = [];
for i = 0:S.NumOfSteps
	eval(['Matrx3D_pso(:, :, i + 1) = S.particle_position_3D_step_', num2str(i), '(:, :);']);
end
N_ = S.NumOfSteps; clear S

for i = 0:1
	title(['DFN flow (mhfem); step NO = ', num2str(i)]);
	p_s = scatter3(Matrx3D_pso(:, 1, i + 1), Matrx3D_pso(:, 2, i + 1), Matrx3D_pso(:, 3, i + 1), 'k', 'o', 'filled');
% 	if (i == 0); pause; else; pause(0.01); end;
% 	if (i ~= N_); delete(p_s); end
end

hold on; figure(1)
Intersection_ = [-0.000000778288324454479152336717, 8.250000000000000000000000000000, -0.000001348034857073798775672913];
scatter3(Intersection_(1), Intersection_(2), Intersection_(3), ...
    'green', 'o', 'filled'); hold on
D1 = [-0.500000000000000000000000000000, 0.000000000000000000000000000000, -0.866025388240814208984375000000];
D2 = [0.642787694931030273437500000000, 0.000000000000000000000000000000, -0.766044378280639648437500000000];

quiver3(Intersection_(1), Intersection_(2), Intersection_(3), ...
    D1(1), D1(2), D1(3), 2, 'linewidth', 1.5); hold on
quiver3(Intersection_(1), Intersection_(2), Intersection_(3), ...
    D2(1), D2(2), D2(3), 2, 'linewidth', 1.5); hold on