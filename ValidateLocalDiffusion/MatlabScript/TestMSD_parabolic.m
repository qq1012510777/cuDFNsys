clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));
deltaT = h5read([currentPath, ...
    '/ParticlePositionResult/DispersionInfo.h5'], "/Delta_T");
Dispersion_local = h5read([currentPath, ...
    '/ParticlePositionResult/DispersionInfo.h5'], "/Dispersion_local");
NumSteps = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], "/NumOfSteps");
%NumSteps = 64
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

f = fittype('a*x^2 + b*x + c', 'independent', 'x', 'coefficients', {'a', 'b', 'c'});
[cfun, goodness] = fit(A', B, f, 'startpoint', [(B(end) - B(1))/(A(end)-A(1)), 0, 0]);


hold on
plot([min_:1: max_], [min_:1: max_] .^2 .* cfun.a + cfun.b .* [min_:1: max_] + cfun.c, 'r-');

ElementAperture = h5read([currentPath, '/MHFEM_1.h5'], '/ElementAperture');
ElementAperture = unique(ElementAperture);
varV = var(ElementAperture .^ 2 / 12 .*20 / 30);
disp(['-----------------'])
disp('fitting mode is: a*x^2 + b*x + c')
disp(['-----------------'])
disp(['goodness.rsquare=', num2str(goodness.rsquare)]);
disp(['a =', num2str(cfun.a)]);
disp(['b =', num2str(cfun.b)]);
disp(['c =', num2str(cfun.c)]);
disp(['-----------------'])
disp(['2 * Dm=', num2str(2 * Dispersion_local)]);
disp(['var[v] = ', num2str(varV)])