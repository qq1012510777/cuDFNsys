clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));

L_d = h5read([currentPath, '/Fractures.h5'], ['/L']);
PercoDir = h5read([currentPath, '/mesh.h5'], ['/group_mesh/Dir']);
DomainDimensionRatio = h5read([currentPath, '/Fractures.h5'], ['/DomainDimensionRatio']);
L = L_d * DomainDimensionRatio(PercoDir + 1);

deltaL = 50;

DeltaT = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], ['/Delta_T']);

D = [];

Sa = [];

QE = [deltaL:deltaL:L];

for i = QE
    data = [];
    if(i ~= L)
        data = h5read([currentPath, '/ParticlePositionResult/ParticlePosition_WhichStepDoesTheParticleReachedControlPlane_', num2str(i, '%.6f'), '_m.h5'], ...
            ['/TravelTime_ReachedControlPlane_', num2str(i, '%.6f'), '_m']);
    else
        data = h5read([currentPath, '/ParticlePositionResult/ParticlePosition_WhichStepDoesTheParticleReached.h5'], ...
            ['/WhichStepDoesTheParticleReached']);
    end

    AS = find(data(:, 1) == -1);

    data(AS, :) = [];
    
    FPT = data(:, 1);
    ParticleTrajectory = data(:, 2);
    
    FPT = double(FPT) .* deltaL;
    
    T1 = mean(FPT);
    
    T2 = sum(FPT.^2) / size(FPT, 1);
    
    % D = [D; (T2 / (T1^2) - 1) * i * (i / T1) / 2];
    D = [D; mean((FPT - T1).^2)/i];
    
    Sa = [Sa; mean(ParticleTrajectory)];
end

figure(1)
scatter([QE], D, 'o')
hold on
xlabel('Linear distence')
ylabel('Quantification of the dispersion')

figure(2)
scatter([QE], Sa./[QE]', 'o')
hold on
xlabel('Linear distence')
ylabel('Particle trajectory (m)')