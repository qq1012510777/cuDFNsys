clc;
close all;
clear all;
currentPath = fileparts(mfilename('fullpath'));
L = h5read([currentPath, '/DFN_VISUAL.h5'], '/L_m');
DomainDimensionRatio = h5read([currentPath, '/DFN_VISUAL.h5'], '/DomainDimensionRatio');
cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, L, -L; -L,L, L];
cube_frame(:, 1) = 0.5 .* cube_frame(:, 1) .* DomainDimensionRatio(1); cube_frame(:, 2) = 0.5 .* cube_frame(:, 2) .* DomainDimensionRatio(2); cube_frame(:, 3) = 0.5 .* cube_frame(:, 3) .* DomainDimensionRatio(3);
figure(1); view(3); title('Discete fracture network'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on
patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on
Frac_NUM_verts = h5read([currentPath, '/DFN_VISUAL.h5'], '/Frac_NUM_verts');
verts = h5read([currentPath, '/DFN_VISUAL.h5'], '/verts');

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

intersections = h5read([currentPath, '/DFN_VISUAL.h5'], '/intersections');
if (size(intersections, 2) == 1); intersections=intersections'; end;
intersections_verts = [intersections(:, [1 2 3]); intersections(:, [4 5 6])];
intersections_structures = [[1:size(intersections, 1)]', [1:size(intersections, 1)]' + size(intersections, 1)];
patch('Vertices', intersections_verts, 'Faces', intersections_structures, 'FaceVertexCData', ones(size(intersections_verts, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0, 'linewidth', 3, 'edgecolor', 'r'); view(3); colorbar; hold on;
axis([-1.1 / 2 * DomainDimensionRatio(1) * L,  1.1 / 2 * DomainDimensionRatio(1) * L, -1.1 / 2 * DomainDimensionRatio(2) * L, 1.1 / 2 * DomainDimensionRatio(2) * L, -1.1 / 2 * DomainDimensionRatio(3) * L, 1.1 / 2 * DomainDimensionRatio(3) * L]);
pbaspect([DomainDimensionRatio]); hold on

NumClusters = h5read([currentPath, '/DFN_VISUAL.h5'], '/NumClusters');
ListClusters = {};
for i=1:NumClusters
	ListClusters{i, 1} = h5read([currentPath, '/DFN_VISUAL.h5'], ['/Cluster_', num2str(i)]);
end
PercolationClusters = h5read([currentPath, '/DFN_VISUAL.h5'], '/PercolationClusters');
figure(2); title('DFN highlighting clusters'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); view(3); hold on
axis([-1.1 / 2 * DomainDimensionRatio(1) * L,  1.1 / 2 * DomainDimensionRatio(1) * L, -1.1 / 2 * DomainDimensionRatio(2) * L, 1.1 / 2 * DomainDimensionRatio(2) * L, -1.1 / 2 * DomainDimensionRatio(3) * L, 1.1 / 2 * DomainDimensionRatio(3) * L]);
pbaspect([DomainDimensionRatio]); hold on
patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on
colorValue = zeros(size(element, 1), 1);
for i = 1:size(ListClusters, 1)
	KK = ListClusters{i, 1};
	if(ismember(PercolationClusters, i))
		colorValue(KK) = 1.2;
	else
		colorValue(KK) = rand;
	end
end
patch('Vertices', verts, 'Faces', element, 'FaceVertexCData', colorValue, 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 1); view(3); colorbar; hold on;
colorbar;

figure(3);
Polar_Orientation = h5read([currentPath, '/DFN_VISUAL.h5'], '/Polar_Orientation');
if (size(Polar_Orientation, 2) == 1 && size(Polar_Orientation, 1) == 2); Polar_Orientation=Polar_Orientation'; end;
polarscatter(Polar_Orientation(:, 1), Polar_Orientation(:, 2), 's', 'filled'); rlim([0 0.5*pi]);
rticks([pi / 12, 2 * pi / 12, 3 * pi / 12, 4 * pi / 12, 5 * pi / 12, 6 * pi / 12 ]);
title(['Fractures', '''',' orientations']); hold on
set(gca,'rticklabel',[]);


%if R values have a lognormal distribution, uncommect the following-------
%% R = h5read([currentPath, '/DFN_VISUAL.h5'], '/R');
%% figure(4);
%% nbins = 20;
%% histfit(R, nbins, 'lognormal');
%% pd=fitdist(R,'lognormal')
%% Ex = exp(pd.mu + pd.sigma^2*0.5)
%% Dx = exp(2*pd.mu+pd.sigma^2)*(exp(pd.sigma^2)-1)
