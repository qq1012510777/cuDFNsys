clc
clear all
close all
load('DFN_II.mat');
L = 0.5 * 30;
cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, L, -L; -L,L, L];
figure(1); view(3); title('Discete fracture network'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on
patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on

MAX_num_fracs = max(Frac_NUM_verts);
NUM_fracs = size(Frac_NUM_verts, 1);
element = NaN(NUM_fracs, MAX_num_fracs);

tmpcc = 1;
for i = 1:NUM_fracs
	tmpkk = Frac_NUM_verts(i);
	for j = 1:tmpkk
		element(i, j) = tmpcc; tmpcc = tmpcc + 1;
	end
end

patch('Vertices', verts, 'Faces', element, 'FaceVertexCData', verts(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 1); view(3); colorbar; hold on;


patch('Vertices', verts, 'Faces', element, 'FaceVertexCData', verts(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 1); view(3); colorbar; hold on;
hold on
A = [-0.000000000000003936225248103212, 14.977400194053657855874917004257, -0.000000000000006817742119750171];

C = [0.024055000000000006932232565759, 14.928589999999999804458639118820, 0.041664482176069342345137158645];

D = [-0.014072818670572912261307330084, 14.912456485088345203848803066649, 0.016771332208570625954546784442];

R = [0.014072818670565033147279443426, 15.042343903018970507901030941866, -0.016771332208584257411621010192];

%--------------

vec_ac = (C - A) .* 200;
C = A + vec_ac; 

vec_ad = (D - A) .* 200;
D = A + vec_ad;

K = [A; C];

F = [A; D];

vec_ar = (R - A) .* 200;
R = A + vec_ar; 
G = [A; R];

figure(1)
view(3)
xlabel('x');
ylabel('y')
zlabel('z'); hold on
plot3(K(:, 1), K(:, 2), K(:, 3), 'r-*', 'linewidth', 2);
hold on
plot3(F(:, 1), F(:, 2), F(:, 3), 'b-*', 'linewidth', 2); hold on
plot3(G(:, 1), G(:, 2), G(:, 3), 'g-*', 'linewidth', 2); hold on