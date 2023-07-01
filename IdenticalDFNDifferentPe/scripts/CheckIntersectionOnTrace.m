clc
clear all
close all

figure(1)
scatter(4.634798498683654344176829908974, -1.944540988114024315791539265774, 'o');
hold on

Edge_ = [[-1.198854374732889249699496758694, 6.263827635840454099991347902687];
          [-1.392475241778327488262334554747, 6.300123122240034234664562973194]];

plot(Edge_(:, 1), Edge_(:, 2), '+-'); hold on
