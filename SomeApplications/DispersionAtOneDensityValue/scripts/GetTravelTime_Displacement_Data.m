clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));

ModelNo = [1:20];

% structYY = {'DeltaT', [], 'FPT', [], 'Displacement', [], ...
%               'L', [], 'NumFractures', []};

for i = ModelNo

    if not(isfolder([currentPath, '/DFN_', num2str(i)]))
        disp(['DFN_', num2str(i), ' is not existing'])
        continue
    end

    data_ = h5read([currentPath, '/DFN_', num2str(i), '/ParticlePositionResult/ParticlePosition_WhichStepDoesTheParticleReached.h5'], ['/WhichStepDoesTheParticleReached']);

    DeltaT = h5read([currentPath, '/DFN_', num2str(i), '/ParticlePositionResult/DispersionInfo.h5'], ['/Delta_T']);

    FPT = data_(:, 1);

    Displacement = data_(:, 2);

    clear data_

    AS = find(FPT == -1);

    FPT(AS) = [];

    Displacement(AS) = [];

    clear AS

%     [Frequency, Data_bin] = hist(FPT, 2000);
% 
%     figure(1)
%     subplot(1, 2, 1)
%     loglog(Data_bin, Frequency, '-'); hold on
% 
% 
%     [Frequency, Data_bin] = hist(Displacement, 400);
%     subplot(1, 2, 2)
%     loglog(Data_bin, Frequency, '-'); hold on

    i
    disp(mean(FPT))
    disp(mean(Displacement))
%     pause()
%     close 1
end
