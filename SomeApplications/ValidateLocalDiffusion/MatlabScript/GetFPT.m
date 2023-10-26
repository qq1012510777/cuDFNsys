clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));

deltaT = h5read([currentPath, ...
    '/DispersionInfo.h5'], "/Delta_T");

Dispersion_local = h5read([currentPath, ...
    '/DispersionInfo.h5'], "/Dispersion_local");

FPT = h5read([currentPath, '/ParticlePosition_WhichStepDoesTheParticleReached.h5'], "/WhichStepDoesTheParticleReached");
FPT = FPT(:, 1);
AS = find(FPT == -1);

FPT(AS) = [];
FPT = double(FPT .* deltaT);

nbines = 120;
    
[Frequency, Data_bin] = hist(FPT, nbines);

figure(1)
F(1) = scatter(Data_bin, Frequency, 'o');hold on


% t = [min(FPT):(max(FPT)-min(FPT))/120:max(FPT)];
% 
% D = 2.22e-8;
% u = 2.22e-7;
% 
% Fx = size(FPT, 1) ./ (4 * pi * D .* t) .* exp(-(30 - t .* u).^2 ./ (4 .* D .* t));
% 
% F(2) = plot(t, Fx, 'r-');
% 
% xlabel('Time (s)'); ylabel('Frequency');
% 
% legend([F], 'Random walk', 'Analytical solution')