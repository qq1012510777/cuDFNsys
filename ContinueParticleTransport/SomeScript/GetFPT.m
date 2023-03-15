clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));

deltaT = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], ...
    ['/Delta_T']);

FPT = h5read([currentPath, '/ParticlePositionResult/ParticlePosition_WhichStepDoesTheParticleReached.h5'], ...
    '/WhichStepDoesTheParticleReached');

AS = find(FPT == -1);

FPT(AS) = [];
FPT = double(FPT) .* deltaT;

nbines = 40;

[Frequency, Data_bin] = hist(FPT, nbines);

figure(1)
loglog(Data_bin, Frequency, '-');

hold on
xlabel('Travel time')
ylabel('Frequency')
