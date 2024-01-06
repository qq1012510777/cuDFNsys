clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));
Filepathss = currentPath;
kklo = find(Filepathss == '/');
Filepathss(kklo(end):end) = [];

Theta_s = [10, 30, 45];

title_matrix = ["(a)", "(b)", "(c)"; "(e)", "(f)", "(g)"];

for i = 1:3
    Filepathss_tt = [Filepathss, '/SingleFrac', char(num2str(i))];

    coordinate_3D = h5read([Filepathss_tt, '/DFN_FLOW_VISUAL.h5'], '/coordinate_3D');
    element_3D = h5read([Filepathss_tt, '/DFN_FLOW_VISUAL.h5'], '/element_3D');
    velocity_center_grid = h5read([Filepathss_tt, '/DFN_FLOW_VISUAL.h5'], '/velocity_center_grid');
    pressure_eles = h5read([Filepathss_tt, '/DFN_FLOW_VISUAL.h5'], '/pressure_eles');
    VelocityNorm = vecnorm(velocity_center_grid');
    
    L = h5read([Filepathss_tt, '/DFN_FLOW_VISUAL.h5'], '/L_m');
    DomainDimensionRatio = h5read([Filepathss_tt, '/DFN_FLOW_VISUAL.h5'], '/DomainDimensionRatio');
    cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L,-L; -L, L, -L; -L,L, L];
    cube_frame(:, 1) = 0.5 .* cube_frame(:, 1) .* DomainDimensionRatio(1); cube_frame(:, 2) = 0.5 .* cube_frame(:, 2) .* DomainDimensionRatio(2); cube_frame(:, 3) = 0.5 .* cube_frame(:, 3) .* DomainDimensionRatio(3);
    figure(1); 
    subplot(2, 3, i)
    view(3); 
    title(title_matrix(1, i)); 
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on
    patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4;5 6 7 8;9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on
    
    axis([-1.1 / 2 * DomainDimensionRatio(1) * L,  1.1 / 2 * DomainDimensionRatio(1) * L, -1.1 / 2 * DomainDimensionRatio(2) * L, 1.1 / 2 * DomainDimensionRatio(2) * L, -1.1 / 2 * DomainDimensionRatio(3) * L, 1.1 / 2 * DomainDimensionRatio(3) * L]);
    pbaspect([DomainDimensionRatio]); hold on
    patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', pressure_eles, 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 1); view(3); hold on
    Cb = colorbar;
    Cb.Title.String = 'Hydraulic head [L]';
    Cb.Location = 'south';
    caxis([0, 30]);
    font_size_t=20;
    set(gca, 'FontSize', font_size_t);
    
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
    
    subplot(2, 3, i+3)
    title(title_matrix(2, i)); hold on
    scatter(A, B, 'o', 'MarkerEdgeColor', 'k'); hold on
    scatter(A, A .* 2 .* Dispersion_local .* (cos(Theta_s(i) / 180 * pi)) .^ 2 + B(1), '^', 'MarkerEdgeColor', 'r'); hold on
    
    set(gca, 'FontSize', font_size_t);
    ylabel('$\sigma_L (t)$', 'FontSize', font_size_t, 'interpreter', 'latex');
    xlabel('$t$', 'FontSize', font_size_t, 'interpreter', 'latex');
end

