clear all 
close all
clc

MSD = [];

StepNUM = 6003;
for i = 0:StepNUM
    HJ = h5read(['Dispersion_MeanSquareDisplacement.h5'], ['/MSD_', num2str(i,'%010d')]);
    MSD = [MSD, HJ];
    
end

figure(1)
title('MSD vs. time')
xlabel('Time step', 'interpreter', 'latex')
ylabel('$\left<z^2 (t)\right>$', 'interpreter', 'latex')
hold on

scatter(0:StepNUM, MSD)

hold on; set(gca, 'FontSize', 16);