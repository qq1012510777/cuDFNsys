clear
clc
close all

Q = qGetRotQuaternion(45 * pi / 180, [0, 0, 1]);

Prot = qRotatePoint( [5 0 0], Q )