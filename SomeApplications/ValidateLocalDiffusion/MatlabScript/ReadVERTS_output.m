clc;
close all;
clear all;
currentPath = fileparts(mfilename('fullpath'));
FG= find(currentPath == '/');
FilePathas_ = currentPath(1:FG(end));
FilePathas_(end) = [];
verts = h5read([FilePathas_, '/DFN_II.h5'], '/verts');

NormalVec = cross(verts([2:4:end], :) - verts([1:4:end], :), ...
    verts([4:4:end], :) - verts([1:4:end], :), 2);
NormalVec = [NormalVec' ./ vecnorm(NormalVec')]';
NormalVec(:, 1) = NormalVec(:, 1) .* sign(NormalVec(:, 3));
NormalVec(:, 2) = NormalVec(:, 2) .* sign(NormalVec(:, 3));
NormalVec(:, 3) = NormalVec(:, 3) .* sign(NormalVec(:, 3));

Centerse = 0.5 .* (verts([3:4:end], :) + verts([1:4:end], :));

String_ss = [repmat(['Fracture_'], size(NormalVec, 1), 1), cell2mat(cellstr(num2str([1:size(NormalVec, 1)]', '%03d')))];

String_ss = [String_ss, repmat([','],size(NormalVec, 1), 1)];

String_ss = [String_ss, cell2mat(cellstr(num2str(NormalVec(:, 1), '%.7f'))), repmat([','],size(NormalVec, 1), 1)];
String_ss = [String_ss, cell2mat(cellstr(num2str(NormalVec(:, 2), '%.7f'))), repmat([','],size(NormalVec, 1), 1)];
String_ss = [String_ss, cell2mat(cellstr(num2str(NormalVec(:, 3), '%.7f'))), repmat([','],size(NormalVec, 1), 1)];

String_ss = [String_ss, cell2mat(cellstr(num2str(Centerse(:, 1), '%.7f'))), repmat([','],size(NormalVec, 1), 1)];
String_ss = [String_ss, cell2mat(cellstr(num2str(Centerse(:, 2), '%.7f'))), repmat([','],size(NormalVec, 1), 1)];
String_ss = [String_ss, cell2mat(cellstr(num2str(Centerse(:, 3), '%.7f'))), repmat([','],size(NormalVec, 1), 1)];

String_ss = [String_ss, repmat(['70,'],size(NormalVec, 1), 1)];

String_ss = [String_ss, cell2mat(cellstr(num2str(0.1 + rand(size(NormalVec, 1), 1) .* 0.4, '%.20f'))), repmat([','],size(NormalVec, 1), 1)];
String_ss = [String_ss, cell2mat(cellstr(num2str(5e-10 + rand(size(NormalVec, 1), 1) .* 4.9995e-06, '%.20f'))), repmat([','],size(NormalVec, 1), 1)];

disp(String_ss)