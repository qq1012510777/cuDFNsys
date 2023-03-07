clc
clear all
close all

deltaT = 5e4;

FPT = h5read("ParticlePosition_WhichStepDoesTheParticleReached.h5", "/WhichStepDoesTheParticleReached");

AS = find(FPT == -1);

FPT(AS) = [];
FPT = double(FPT .* deltaT);

nbines = 400;
    
[Frequency, Data_bin] = hist(FPT, nbines);

figure(1)
loglog(Data_bin(1:100), Frequency(1:100), '-');