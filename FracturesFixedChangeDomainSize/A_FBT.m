clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));

N_steps = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], ['/NumOfSteps']);
N_particles = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], ['/NumParticles']);
delta_t = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], ['/Delta_T']);
BlockNOPresent = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], '/BlockNOPresent');
SizeOfDataBlock = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], '/SizeOfDataBlock');

BT_time = zeros(N_particles, 1) - 1;

clear S;

for i = [N_steps:-1:1]
    H5name = [currentPath, '/ParticlePositionResult/ParticlePositionBlock', num2str(ceil(double(i) / double(SizeOfDataBlock)), '%010d'), '.h5'];
	S = h5read(H5name, ['/ParticleIDAndElementTag_', num2str(i, '%010d')]);
	ExistingParticle = S(:, 1) + 1.0;
    AS = [1:1:N_particles];
    AS(ismember(AS, ExistingParticle)) = [];
    BT_time(AS, :) = double(i) * delta_t;
    
    clear S AS ExistingParticle
    disp([num2str(i)]);
end

BT_time2 = find(BT_time ~= -1);

BT_time = BT_time(BT_time2, :);

% histogram(BT_time)

[ycdf,xcdf] = cdfcalc(BT_time);
% plot(xcdf, ycdf(1:end-1)); hold on

yccdf = 1-ycdf(1:end-1);
figure(1)
title('CCDF of particle trval time', 'interpreter', 'latex')
loglog(xcdf,yccdf, 'o-r'); hold on
ylabel('CCDF', 'interpreter', 'latex');
xlabel('$\Delta t$', 'interpreter', 'latex')

nbines = 30;
figure(2)
[Frequency, Data_bin]=hist(BT_time, nbines);
loglog(Data_bin,Frequency, 'ok'); hold on
title('Distribution of particle trval time', 'interpreter', 'latex')
ylabel('Frequency', 'interpreter', 'latex');
xlabel('$\Delta t$', 'interpreter', 'latex')
