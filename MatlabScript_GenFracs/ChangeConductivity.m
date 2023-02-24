clc
clear all
close all
currentPath = fileparts(mfilename('fullpath'));
NumFracs = h5read([currentPath, '/Fractures.h5'], '/NumFractures');
L = h5read([currentPath, '/Fractures.h5'], '/L');

gamma = 5e-4;
beta = 0.1;

new_filename = 'Fractures_new.h5';

h5create([currentPath, '/', new_filename], '/L', size(L));
h5write([currentPath, '/', new_filename], '/L', L)

h5create([currentPath, '/', new_filename], '/NumFractures', size(NumFracs));
h5write([currentPath, '/', new_filename], '/NumFractures', NumFracs);

for i = 1:NumFracs
    Center = h5read([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/Center']);
    ConnectModelSurf = h5read([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/ConnectModelSurf']);
    NormalVec = h5read([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/NormalVec']);
    NumVertsTruncated = h5read([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/NumVertsTruncated']);
    Radius = h5read([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/Radius']);
    Verts3D = h5read([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/Verts3D']);
    Verts3DTruncated = h5read([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/Verts3DTruncated']);

    Conductivity = 1.0e-8; %gamma^3.0 / 12 * Radius^(3.0*beta);

    h5create([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/Conductivity'], size(Conductivity));
    h5create([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/Verts3D'], size(Verts3D));
    h5create([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/Center'], size(Center));
    h5create([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/Verts3DTruncated'], size(Verts3DTruncated));
    h5create([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/NumVertsTruncated'], size(NumVertsTruncated));
    h5create([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/Radius'], size(Radius));
    h5create([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/ConnectModelSurf'], size(ConnectModelSurf));
    h5create([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/NormalVec'], size(NormalVec));

    h5write([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/Conductivity'], (Conductivity));
    h5write([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/Verts3D'], (Verts3D));
    h5write([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/Center'], (Center));
    h5write([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/Verts3DTruncated'], (Verts3DTruncated));
    h5write([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/NumVertsTruncated'], (NumVertsTruncated));
    h5write([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/Radius'], (Radius));
    h5write([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/ConnectModelSurf'], (ConnectModelSurf));
    h5write([currentPath, '/', new_filename], ['/Fracture_', num2str(i), '/NormalVec'], (NormalVec));
end
