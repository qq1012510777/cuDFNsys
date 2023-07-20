clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));

deltaT = h5read([currentPath, ...
    '/ParticlePositionResult/DispersionInfo.h5'], "/Delta_T");

Dispersion_local = h5read([currentPath, ...
    '/ParticlePositionResult/DispersionInfo.h5'], "/Dispersion_local");

FPT = h5read([currentPath, '/ParticlePositionResult/ParticlePosition_WhichStepDoesTheParticleReached.h5'], "/WhichStepDoesTheParticleReached");
FPT = FPT(:, 1);
AS = find(FPT == -1);

FPT(AS) = [];
%FPT = double(FPT .* deltaT);
FPT = double(FPT);
[F, X] = ecdf(FPT);

figure(1)
scatter(X, F, 'o');hold on

%  exp(V .* 30 ./ D ) .* erfc((30 + V.*t)./(2 .* sqrt(D .* t))) 
f = fittype('0.5 .* ( erfc((30 - V .* t)./(2 .* sqrt(D .* t))) +exp(V .* 30 ./ D ) .* erfc((30 + V.*t)./(2 .* sqrt(D .* t)))  )', 'independent', 't', 'coefficients', {'V', 'D'});
[cfun, goodness] = fit(X, F, f, 'startpoint', [0.2, 0.01]);

V = cfun.V;
D = cfun.D;

t = X;
F2 = 0.5 .* ( erfc((30 - V .* t)./(2 .* sqrt(D .* t))) +exp(V .* 30 ./ D ) .* erfc((30 + V.*t)./(2 .* sqrt(D .* t)))  );

figure(1)
plot(X, F2, '-');
disp(['V_t = ', num2str(V), ', D_t = ', num2str(D), ', goodness = ', num2str(goodness.rsquare)])
% nbines = 25;

%  D = 8e-9;
%  V= 2e-10;
%  Tstart = 6e9;
%  Tfinish=1e12;
%  t = [Tstart:((Tfinish-Tstart) / 100): Tfinish];
%  % exp(V .* 30 ./ D ) .* erfc((30 + V.*t)./(2 .* sqrt(D .* t))) 
%  G = 0.5 .* ( erfc((30 - V .* t)./(2 .* sqrt(D .* t))) +0);
%  figure(2)
%  scatter(t, G)
%     
% [Frequency, Data_bin] = hist(FPT, nbines);
% %Frequency = Frequency ./sum(Frequency);
% 
% figure(1)
% title("FPT curve")
% xlabel("time")
% ylabel("frequency")
% hold on
% F(1) = scatter(Data_bin, Frequency, 'o');hold on
% 
% f = fittype('1 / (4 * pi * D * t) * exp(-(28 - t * V)^2/(4*D*t))', 'independent', 't', 'coefficients', {'V', 'D'});
% [cfun, goodness] = fit(Data_bin', Frequency', f, 'startpoint', [1e-9, 2e-9]);
% 
% goodness.rsquare
