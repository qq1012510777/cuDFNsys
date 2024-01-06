clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));
Filepathss = currentPath;
kklo = find(Filepathss == '/');
Filepathss(kklo(end):end) = [];

Filepathss_tt = [Filepathss, '/NumerousInclinedFractures'];
font_size_t = 16;

element_Frac_Tag = h5read([Filepathss_tt, '/DFN_MESH_VISUAL.h5'], '/element_Frac_Tag');
num_frac = max(element_Frac_Tag);

velocity_center_grid = h5read([Filepathss_tt, '/DFN_FLOW_VISUAL.h5'], '/velocity_center_grid');
VelocityNorm = vecnorm(velocity_center_grid');
ElementAperture = h5read([Filepathss_tt, '/DFN_FLOW_VISUAL.h5'], '/ElementAperture');
VelocityNorm = [vecnorm(velocity_center_grid')]' ./ ElementAperture;
pressure_eles = h5read([Filepathss_tt, '/DFN_FLOW_VISUAL.h5'], '/pressure_eles');

SizeOfDataBlock = h5read([Filepathss_tt, '/ParticlePositionResult/DispersionInfo.h5'], '/SizeOfDataBlock');

coordinate_3D = h5read([Filepathss_tt, '/DFN_FLOW_VISUAL.h5'], '/coordinate_3D');
element_3D = h5read([Filepathss_tt, '/DFN_FLOW_VISUAL.h5'], '/element_3D');
L = h5read([Filepathss_tt, '/DFN_FLOW_VISUAL.h5'], '/L_m');
DomainDimensionRatio = h5read([Filepathss_tt, '/DFN_FLOW_VISUAL.h5'], '/DomainDimensionRatio');
cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, L, -L; -L,L, L];
cube_frame(:, 1) = 0.5 .* cube_frame(:, 1) .* DomainDimensionRatio(1); cube_frame(:, 2) = 0.5 .* cube_frame(:, 2) .* DomainDimensionRatio(2); cube_frame(:, 3) = 0.5 .* cube_frame(:, 3) .* DomainDimensionRatio(3);
figure(1); 
subplot(2, 2, 1);
view(3); 
title('(a)'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on
patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on

axis([-1.1 / 2 * DomainDimensionRatio(1) * L,  1.1 / 2 * DomainDimensionRatio(1) * L, -1.1 / 2 * DomainDimensionRatio(2) * L, 1.1 / 2 * DomainDimensionRatio(2) * L, -1.1 / 2 * DomainDimensionRatio(3) * L, 1.1 / 2 * DomainDimensionRatio(3) * L]);
pbaspect([DomainDimensionRatio]); hold on
patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', pressure_eles, 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 1); view(3); hold on
Cb = colorbar;
% Cb.Title.String = 'Hydraulic head [L]';
caxis([0, 30]);
set(gca, 'FontSize', font_size_t);

figure(1)
subplot(2, 2, 2);
title('(b)'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on
patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on
view(3)
for i = 0
	H5name = []; H5name_2D = [];
	if (i == 0); H5name = [Filepathss_tt, '/ParticlePositionResult/ParticlePositionInit_3D.h5']; else; H5name = [Filepathss_tt, '/ParticlePositionResult/ParticlePositionBlock', num2str(ceil(double(i) / double(SizeOfDataBlock)), '%010d'), '_3D.h5']; end;
	if (i == 0); H5name_2D = [Filepathss_tt, '/ParticlePositionResult/ParticlePositionInit.h5']; else; H5name_2D = [Filepathss_tt, '/ParticlePositionResult/ParticlePositionBlock', num2str(ceil(double(i) / double(SizeOfDataBlock)), '%010d'), '.h5']; end;
	S = h5read(H5name, ['/Step_', num2str(i, '%010d')]);
	ParticleID = h5read(H5name_2D, ['/ParticleIDAndElementTag_', num2str(i, '%010d')]);
	Matrx3D_pso = [];%NaN(N_particles, 3);
	Matrx3D_pso(:, :) = S(:, [1 2 3]);
	p_s = scatter3(Matrx3D_pso(:, 1), Matrx3D_pso(:, 2), Matrx3D_pso(:, 3), [], ParticleID(:, 1), 'filled'); clear Matrx3D_pso

end
set(gca, 'FontSize', font_size_t);

figure(1)
subplot(2, 2, 3);
title('(c)'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on
patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on
view(3)
for i = 151
	H5name = []; H5name_2D = [];
	if (i == 0); H5name = [Filepathss_tt, '/ParticlePositionResult/ParticlePositionInit_3D.h5']; else; H5name = [Filepathss_tt, '/ParticlePositionResult/ParticlePositionBlock', num2str(ceil(double(i) / double(SizeOfDataBlock)), '%010d'), '_3D.h5']; end;
	if (i == 0); H5name_2D = [Filepathss_tt, '/ParticlePositionResult/ParticlePositionInit.h5']; else; H5name_2D = [Filepathss_tt, '/ParticlePositionResult/ParticlePositionBlock', num2str(ceil(double(i) / double(SizeOfDataBlock)), '%010d'), '.h5']; end;
	S = h5read(H5name, ['/Step_', num2str(i, '%010d')]);
	ParticleID = h5read(H5name_2D, ['/ParticleIDAndElementTag_', num2str(i, '%010d')]);
	Matrx3D_pso = [];%NaN(N_particles, 3);
	Matrx3D_pso(:, :) = S(:, [1 2 3]);
	p_s = scatter3(Matrx3D_pso(:, 1), Matrx3D_pso(:, 2), Matrx3D_pso(:, 3), [], ParticleID(:, 1), 'filled'); clear Matrx3D_pso

end
set(gca, 'FontSize', font_size_t);

figure(1)
subplot(2, 2, 4);
deltaT = h5read([Filepathss_tt, ...
                 '/ParticlePositionResult/DispersionInfo.h5'], "/Delta_T");
Dispersion_local = h5read([Filepathss_tt, ...
                           '/ParticlePositionResult/DispersionInfo.h5'], "/Dispersion_local");
NumSteps = h5read([Filepathss_tt, '/ParticlePositionResult/DispersionInfo.h5'], "/NumOfSteps");
% NumSteps = 30
Mean_cMSD_MSD = [];

for j = 0:NumSteps

    try
        AE = h5read([Filepathss_tt, '/Dispersion_MeanSquareDisplacement.h5'], ['/Variance_of_xyz_', num2str(j, '%010d')]);
        Mean_cMSD_MSD = [Mean_cMSD_MSD; AE'];
    catch
        NumSteps = size(Mean_cMSD_MSD, 1);
        break
    end

end

A = double(deltaT) .* double([0:NumSteps]);
B = [Mean_cMSD_MSD(:, 3)];

sq = [];
sq(1) = scatter(A, B, 'o', 'MarkerEdgeColor', 'k'); hold on
% scatter(A, A .* 2 .* Dispersion_local .* (cos(Theta_s(i) / 180 * pi)) .^ 2 + B(1), '^', 'MarkerEdgeColor', 'r'); hold on
set(gca, 'FontSize', font_size_t);
ylabel('$\sigma_L (t)$', 'FontSize', font_size_t, 'interpreter', 'latex');
xlabel('$t$', 'FontSize', font_size_t, 'interpreter', 'latex');

verts = h5read([Filepathss_tt, '/DFN_VISUAL.h5'], '/verts');

NormalVec = cross(verts([2:4:end], :) - verts([1:4:end], :), ...
    verts([4:4:end], :) - verts([1:4:end], :), 2);
NormalVec = [NormalVec' ./ vecnorm(NormalVec')]';
NormalVec(:, 1) = NormalVec(:, 1) .* sign(NormalVec(:, 3));
NormalVec(:, 2) = NormalVec(:, 2) .* sign(NormalVec(:, 3));
NormalVec(:, 3) = NormalVec(:, 3) .* sign(NormalVec(:, 3));
Theta_s = pi/2.0 - acos(NormalVec(:, 3));


FracVelocity = zeros(num_frac, 1);

for i = 1:num_frac
    ind_frac = find(element_Frac_Tag == i);
    FracVelocity(i) = mean(VelocityNorm(ind_frac));
end

sq(2) = plot(A, A .* 2 .* Dispersion_local .* mean(cos(Theta_s) .^ 2) + B(1) + A .^2 .* var(FracVelocity .* cos(Theta_s)), 'r-', 'LineWidth', 2); hold on

legend([sq], {'Measured in (a)', 'Predicted'})