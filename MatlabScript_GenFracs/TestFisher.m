clc
clear all
close all

directoryName = "./";

cd(directoryName);

figure(1)
Polar_Orientation = h5read(['./DFN_VISUAL.h5'], '/Polar_Orientation');
if (size(Polar_Orientation, 2) == 1 && size(Polar_Orientation, 1) == 2); Polar_Orientation=Polar_Orientation'; end;
polarscatter(Polar_Orientation(:, 1), Polar_Orientation(:, 2), 's', 'filled'); rlim([0 0.5*pi]);
rticks([pi / 12, 2 * pi / 12, 3 * pi / 12, 4 * pi / 12, 5 * pi / 12, 6 * pi / 12 ]);
title(['Fractures', '''',' orientations']); hold on
set(gca,'rticklabel',[]);

phi__ = Polar_Orientation(:, 1);
theta__ = Polar_Orientation(:, 2);

figure(2)
subplot(1, 2, 1)
[Frequency0, edges] = histcounts(phi__, 'NumBins', 130);
Data_bin0 = (edges(1:end-1) + edges(2:end)) / 2;
DeltaX = Data_bin0(2) - Data_bin0(1);
Frequency0 = Frequency0./(sum(Frequency0) * DeltaX);
P(1) = plot(Data_bin0, Frequency0, 'o'); hold on

xlabel("$\phi$ value", 'Interpreter','latex')
ylabel("Frequency")
ylim([0,1])

subplot(1, 2, 2)
[Frequency0, edges] = histcounts(theta__, 'NumBins', 130);
Data_bin0 = (edges(1:end-1) + edges(2:end)) / 2;
DeltaX = Data_bin0(2) - Data_bin0(1);
Frequency0 = Frequency0./(sum(Frequency0) * DeltaX);
P(2) = plot(Data_bin0, Frequency0, 'o'); hold on

xlabel("$\theta$ value", 'Interpreter','latex')
ylabel("Frequency")
ylim([0,4])

kappa = 20;                      % concentration parameter
theta = linspace(0, 0.5*pi, 1000);   % theta from 0 to pi

% Compute PDF of theta for Fisher distribution
% f(θ) = (κ * sinθ / (2 * sinhκ)) * exp(κ cosθ)
pdf_theta = (kappa .* sin(theta)) ./ (2 * sinh(kappa)) .* exp(kappa .* cos(theta));
subplot(1, 2, 2)
plot(theta, pdf_theta, 'LineWidth', 2); hold on