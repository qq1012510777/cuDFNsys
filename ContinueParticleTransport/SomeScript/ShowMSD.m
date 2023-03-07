clear all 
close all
clc


deltaT = 5e4;
StepNUM = 120003;
StepNUM = 13749;
MSD = zeros(1, StepNUM+1);

for i = 0:StepNUM
    HJ = h5read(['Dispersion_MeanSquareDisplacement.h5'], ['/MSD_', num2str(i,'%010d')]);
    MSD(i+1) = HJ;
end

figure(1)

loglog([0:StepNUM] .* deltaT, MSD);
hold on

title('cMSD vs. time')
xlabel('Time step', 'interpreter', 'latex')
ylabel('$\left<\left|z-\left<z\right>\right|^2\right>$', 'interpreter', 'latex')
hold on

hold on; set(gca, 'FontSize', 16);