clc
clear all
close all

Theta_s = 10; % the angle between the fracture and the mean flow direction

currentPath = fileparts(mfilename('fullpath'));
Filepathss = currentPath;
kklo = find(Filepathss == '/');
Filepathss(kklo(end):end) = [];

deltaT = h5read([Filepathss, ...
    '/ParticlePositionResult/DispersionInfo.h5'], "/Delta_T");
Dispersion_local = h5read([Filepathss, ...
    '/ParticlePositionResult/DispersionInfo.h5'], "/Dispersion_local");
NumSteps = h5read([Filepathss, '/ParticlePositionResult/DispersionInfo.h5'], "/NumOfSteps");
% NumSteps = 30
Mean_cMSD_MSD = [];
for i = 0:NumSteps
    try
        AE = h5read([Filepathss, '/Dispersion_MeanSquareDisplacement.h5'], ['/Variance_of_xyz_', num2str(i, '%010d')]);
        Mean_cMSD_MSD = [Mean_cMSD_MSD; AE'];
    catch
        NumSteps = size(Mean_cMSD_MSD, 1);
        break
    end
end
figure(1)
A = double(deltaT) .* double([0:NumSteps]);
B = [Mean_cMSD_MSD(:, 3)];
scatter(A, B)

min_ = min(A');
max_ = max(A');

f = fittype('a*x + c', 'independent', 'x', 'coefficients', {'a', 'c'});
[cfun, goodness] = fit(A', B, f, 'startpoint', [(B(end) - B(1))/(A(end)-A(1)), 0]);


hold on
plot([min_, max_], [min_, max_] .* cfun.a + cfun.c, 'r-');

disp(['goodness.rsquare=', num2str(goodness.rsquare)]);
disp(['slope =', num2str(cfun.a)]);
disp(['Dm_input =', num2str(Dispersion_local)]);
disp(['slope / (Dm_input * cos(theta)^2) = ', num2str(cfun.a/(Dispersion_local .* (cos(Theta_s * pi / 180) .^ 2)))])

figure(2)
scatter(A, B, 'o', 'MarkerEdgeColor', 'k'); hold on
scatter(A, A .* 2 .* Dispersion_local .* (cos(Theta_s / 180 * pi)) .^ 2 + B(1), '^', 'MarkerEdgeColor', 'r'); hold on

font_size_t=20;
set(gca, 'FontSize', font_size_t);
ylabel('$\sigma_L (t)$', 'FontSize', font_size_t, 'interpreter', 'latex');
xlabel('$t$', 'FontSize', font_size_t, 'interpreter', 'latex');
% set(gca, 'XScale', 'log', 'YScale', 'log')