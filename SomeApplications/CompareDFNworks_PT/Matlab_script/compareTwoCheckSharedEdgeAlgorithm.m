clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));
klj = find(currentPath == '/');
FileDir_c = currentPath(1:klj(end) - 1);

A1 = h5read([FileDir_c, '/EdgesSharedEle.h5'], '/data');
A2 = h5read([FileDir_c, '/EdgesSharedEle_III.h5'], '/data');

as = sum(sum(A1-A2));