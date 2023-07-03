clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));
deltaT = h5read([currentPath, ...
    '/ParticlePositionResult/DispersionInfo.h5'], "/Delta_T");
Dispersion_local = h5read([currentPath, ...
    '/ParticlePositionResult/DispersionInfo.h5'], "/Dispersion_local");
NumSteps = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], "/NumOfSteps");

Mean_cMSD_MSD = [];
for i = 1:NumSteps
    try
        AE = h5read([currentPath, '/Dispersion_MeanSquareDisplacement.h5'], ['/Mean_MSD_cMSD_', num2str(i, '%010d')]);
        Mean_cMSD_MSD = [Mean_cMSD_MSD; AE'];
    catch
        NumSteps = size(Mean_cMSD_MSD, 1);
        break
    end
end
figure(1)
A = double(deltaT) .* double([1:NumSteps]);
B = [Mean_cMSD_MSD(:, 3) - Mean_cMSD_MSD(:, 1) .^ 2];
scatter(A, B)

min_ = min(A');
max_ = max(A');

f = fittype('a*x + c', 'independent', 'x', 'coefficients', {'a', 'c'});
[cfun, goodness] = fit(A', B, f, 'startpoint', [(B(end) - B(1))/(A(end)-A(1)), 0]);


hold on
plot([min_, max_], [min_, max_] .* cfun.a + cfun.c, 'r-');

disp(['goodness.rsquare=', num2str(goodness.rsquare)]);
disp(['cfun.a=', num2str(cfun.a)]);
disp(['Dispersion_local=', num2str(Dispersion_local)]);