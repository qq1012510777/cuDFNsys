clear all 
close all
clc

IfPredict = 0;

currentPath = fileparts(mfilename('fullpath'));

deltaT = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], ['/Delta_T']);
StepNUM = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], ['/NumOfSteps']);
StepNUM = double(StepNUM);

MSD = zeros(1, StepNUM+1);

BlockNOPresent = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], '/BlockNOPresent');
SizeOfDataBlock = h5read([currentPath, '/ParticlePositionResult/DispersionInfo.h5'], '/SizeOfDataBlock');

for i = 0:StepNUM
    
    try
        HJ = h5read(['Dispersion_MeanSquareDisplacement.h5'], ['/Variance_of_xyz_', num2str(i,'%010d')]);
    catch
        break
    end
    MSD(i+1) = HJ(3);% - HJ(1).^2;

    % try
    %     H5name = []; H5name_2D = [];
	% 	if (i == 0); H5name = [currentPath, '/ParticlePositionResult/ParticlePositionInit_3D.h5']; else; H5name = [currentPath, '/ParticlePositionResult/ParticlePositionBlock', num2str(ceil(double(i) / double(SizeOfDataBlock)), '%010d'), '_3D.h5']; end;
	% 	if (i == 0); H5name_2D = [currentPath, '/ParticlePositionResult/ParticlePositionInit.h5']; else; H5name_2D = [currentPath, '/ParticlePositionResult/ParticlePositionBlock', num2str(ceil(double(i) / double(SizeOfDataBlock)), '%010d'), '.h5']; end;
    % 
    %     Pos = h5read(H5name, ['/Step_', num2str(i, '%010d')]);
    %     Factor_pe = h5read(H5name_2D, ['/FactorPeriodic', num2str(i, '%010d')]);
    %     Pos(:, 3) = Pos(:, 3) + double(Factor_pe).* 60;
    % 
    %     if (~isempty(find(Factor_pe > 7)))
    %         error("11")
    %     end
    % 
    %     MSD(i+1) = var(Pos(:, 3)); 
    % 
    % catch
    %     break
    % end
end
MSD(i+1:end) = [];


figure(1)

% scatter([0:StepNUM] .* deltaT, cMSD);
StartNO = 1;%ceil(StepNUM / 2);
A = [StartNO:size(MSD, 2)-1] .* deltaT;

B = MSD(StartNO + 1:end);
scatter(A, B);
hold on

title('Variance vs. time')
xlabel('Time step', 'interpreter', 'latex')
ylabel('Variance in the z axis', 'interpreter', 'latex')
hold on

if (IfPredict == 0)
    syms x
    f = fittype('a*x + c', 'independent', 'x', 'coefficients', {'a', 'c'});
    [cfun, goodness] = fit(A', B', f); 
    plot(A, A.*cfun.a + cfun.c, 'r-')
    disp(['goodness.rsquare=', num2str(goodness.rsquare)]);
    disp(['cfun.a=', num2str(cfun.a)]);
    disp(['cfun.c=', num2str(cfun.c)]);
else
    A_p = A;
    B_p =0.238362 + A_p.* 2 * 8e-5 .* mean(cos(1/3*pi)^2);
    plot(A_p, B_p, 'r-'); hold on
end