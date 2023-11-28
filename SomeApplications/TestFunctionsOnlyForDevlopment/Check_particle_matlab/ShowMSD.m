clear all 
close all
clc

currentPath = fileparts(mfilename('fullpath'));

deltaT = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], ['/Delta_T']);
StepNUM = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], ['/NumOfSteps']);
StepNUM = double(StepNUM);

MSD = zeros(3, StepNUM+1);

for i = 0:StepNUM
    
    try
        HJ = h5read(['Dispersion_MeanSquareDisplacement.h5'], ['/Variance_of_xyz_', num2str(i,'%010d')]);
    catch
        break
    end
    MSD(i+1) = HJ(3);% - HJ(1).^2;
end
MSD(i+1:end) = [];

figure(1)

% scatter([0:StepNUM] .* deltaT, cMSD);
A = [0:size(MSD, 2)-1] .* deltaT;
scatter(A, MSD);
hold on

title('Variance vs. time')
xlabel('Time step', 'interpreter', 'latex')
ylabel('Variance in the z axis', 'interpreter', 'latex')
hold on

syms x

f = fittype('a*x + c', 'independent', 'x', 'coefficients', {'a', 'c'});
[cfun, goodness] = fit(A', MSD', f);
disp(['goodness.rsquare=', num2str(goodness.rsquare)]);
disp(['cfun.a=', num2str(cfun.a)]);
disp(['cfun.c=', num2str(cfun.a)]);