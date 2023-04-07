clc
clear all
close all

figure(1)
pbaspect([1, 1, 1]); hold on
plot([0, 100, 100, 0, 0], [0, 0, 100, 100, 0], 'k-', ...
    'linewidth', 1.5); hold on
grid on

xticks([0:10:100])

set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
hold on

for i=[10:10:90]
    plot([i, i], [0, 100], 'k-', 'linewidth', 1.5); hold on
end

for i=[10:10:90]
    plot([0, 100], [i, i], 'k-', 'linewidth', 1.5); hold on
end

NumSites=70;

X = randi([1,10], [NumSites, 1]) .* 10 - 5;
Y = randi([1,10], [NumSites, 1]) .* 10 - 5;

for i = [1:NumSites]
    scatter(X, Y, 200, 'o', 'filled', 'markerfacecolor', 'k'); hold on
end