clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));
deltaT = h5read([currentPath, ...
    '/ParticlePositionResult/DispersionInfo.h5'], "/Delta_T");
Dispersion_local = h5read([currentPath, ...
    '/ParticlePositionResult/DispersionInfo.h5'], "/Dispersion_local");
NumSteps = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], "/NumOfSteps");
NumSteps = 100;
Variance = [];
for i = 1:NumSteps
    try
        AE = h5read([currentPath, '/Dispersion_MeanSquareDisplacement.h5'], ['/Variance_of_xyz_', num2str(i, '%010d')]);
        Variance = [Variance; AE'];
    catch
        NumSteps = size(Variance, 1);
        break
    end
end
figure(1)
A = double(deltaT) .* double([1:NumSteps]);
B = [Variance(:, 3)];
scatter(A, B)

min_ = min(A');
max_ = max(A');

%------------------------variance-------------------------------
Frac_NUM_verts = h5read([currentPath, '/DFN_II.h5'], '/Frac_NUM_verts');
verts = h5read([currentPath, '/DFN_II.h5'], '/verts');

MAX_num_fracs = max(Frac_NUM_verts);
NUM_fracs = size(Frac_NUM_verts, 1);
element = NaN(NUM_fracs, MAX_num_fracs);
L = h5read([currentPath, '/DFN_II.h5'], '/L_m');

tmpcc = 1;
for i = 1:NUM_fracs
	tmpkk = Frac_NUM_verts(i);
	for j = 1:tmpkk
		element(i, j) = tmpcc; tmpcc = tmpcc + 1;
	end
end
cos_theta = zeros(NUM_fracs, 1);

LengthFrac = cos_theta;
for i = 1:NUM_fracs
   L_e = norm(verts(element(i, 2), :) -  verts(element(i, 1), :));
   cos_theta(i, 1) = L / L_e;
   LengthFrac(i, 1) = L_e;
end
%----------------------

f = fittype('a*x + c', 'independent', 'x', 'coefficients', {'a', 'c'});
[cfun, goodness] = fit(A', B, f, 'startpoint', [(B(end) - B(1))/(A(end)-A(1)), 0]);


hold on
plot([min_, max_], [min_, max_] .* cfun.a + cfun.c, 'r-');

disp(['goodness.rsquare=', num2str(goodness.rsquare)]);
disp(['slope =', num2str(cfun.a)]);
disp(['Dm=', num2str(Dispersion_local)]);
disp(['slope / Dm = ', num2str(cfun.a/Dispersion_local)])
disp(['slope / (cos ^2 \theta * Dm) = ', num2str(cfun.a/(cos_theta(1).^2 * Dispersion_local))])