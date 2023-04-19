clear all 
close all
clc

currentPath = fileparts(mfilename('fullpath'));

deltaT = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], ['/Delta_T']);
StepNUM = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], ['/NumOfSteps']);
StepNUM = double(StepNUM);

cMSD = zeros(1, StepNUM+1);

for i = 0:StepNUM
    
    try
        HJ = h5read(['Dispersion_MeanSquareDisplacement.h5'], ['/Mean_MSD_cMSD_', num2str(i,'%010d')]);
    catch
        break
    end
    cMSD(i+1) = HJ(3);% - HJ(1).^2;
end
cMSD(i+1:end) = [];

figure(1)

% scatter([0:StepNUM] .* deltaT, cMSD);
scatter([0:size(cMSD, 2)-1], cMSD);
hold on

title('cMSD vs. time')
xlabel('Time step', 'interpreter', 'latex')
ylabel('$\left<\left|z-\left<z\right>\right|^2\right>$', 'interpreter', 'latex')
hold on

hold on; set(gca, 'FontSize', 16);