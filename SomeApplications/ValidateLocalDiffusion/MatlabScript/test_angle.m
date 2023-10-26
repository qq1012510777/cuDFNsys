clc
clear all
close all

v = [-1, 1];
w = [1, -0.75];
a = atan2(w(2) * v(1) - w(1) * v(2),w(1) * v(1)+w(2)*v(2)) * 180.0 / pi;

b = acos(dot(v, w) / (norm(v)*norm(w))) * 180 / pi;