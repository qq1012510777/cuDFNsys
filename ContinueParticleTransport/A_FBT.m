clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));

S = load([currentPath, '/ParticlePositionResult/ParticlePosition_step_0000000.mat']);
N_steps = S.NumOfSteps;
N_particles = size(S.particle_position_3D_step_0, 1);

BT_time = zeros(N_particles, 1) - 1;

clear S;

for i = [N_steps:-1:1]
    S = load([currentPath, '/ParticlePositionResult/ParticlePosition_step_', num2str(i, '%07d'), '.mat']);
    eval(['AS = find(S.particle_position_3D_step_', num2str(i), '(:, 4) == 1);']);
    BT_time(AS, :) = i * 0.01;
    clear S AS
    i
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

nbines = 50;
figure(2)
[Frequency, Data_bin]=hist(BT_time, nbines);
loglog(Data_bin,Frequency, 'ok'); hold on
title('Distribution of particle trval time', 'interpreter', 'latex')
ylabel('Frequency', 'interpreter', 'latex');
xlabel('$\Delta t$', 'interpreter', 'latex')
