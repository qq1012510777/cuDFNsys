clc
clear all
close all

b = 1e-2; % m

mu = 1.002; % Pa * S
k = 2.5e-13; % m * m
L1 = 0.9; % m
L3 = L1;
L4 = L1;
L23 = 0.2; % m
L24 = L23;

X = 1.8; % m
Y = 1.0;
Z = 0.65;

A = Y * b; % m * m

R1 = mu * L1 / k / A; % Pa * S / m ^ 3
R3 = mu * L3 / k / A;
R4 = mu * L4 / k / A;

R23 = mu * L23 / k / A;
R24 = mu * L24 / k / A;

R_t = R1 + (R24 +R4) * (R23 + R3) / (R24 + R4 + R23 + R3);

Delta_P = 1e6; % Pa

V = X * Y * Z; 

n = A * (L1 + L3 + L4 + L23 + L24) / V

Q = Delta_P / R_t % m ^ 3 / S

q = Q / (Z * Y) % m ^ 2 / S

v = q / n

Traval_time = X / v