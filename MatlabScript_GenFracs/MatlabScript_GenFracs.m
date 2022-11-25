clc
clear all
close all
currentPath = fileparts(mfilename('fullpath'));
addpath(genpath([currentPath, '/Quaternion']));

Frac(1) = struct('Conductivity', [], ...
    'Verts3D', [], ...
    'Center', [], ...
    'Verts3DTruncated', [], ...
    'NumVertsTruncated', [], ...
    'Radius', [], ...
    'ConnectModelSurf', [], ...
    'NormalVec', []);

Frac(1).Conductivity = h5read([currentPath, '/FracturesII.h5'], '/Fracture_99/Conductivity');
Frac(1).Verts3D = h5read([currentPath, '/FracturesII.h5'], '/Fracture_99/Verts3D');
Frac(1).Center = h5read([currentPath, '/FracturesII.h5'], '/Fracture_99/Center');
Frac(1).Verts3DTruncated = h5read([currentPath, '/FracturesII.h5'], '/Fracture_99/Verts3DTruncated');
Frac(1).NumVertsTruncated = h5read([currentPath, '/FracturesII.h5'], '/Fracture_99/NumVertsTruncated');
Frac(1).Radius = h5read([currentPath, '/FracturesII.h5'], '/Fracture_99/Radius');
Frac(1).ConnectModelSurf = h5read([currentPath, '/FracturesII.h5'], '/Fracture_99/ConnectModelSurf');
Frac(1).NormalVec = h5read([currentPath, '/FracturesII.h5'], '/Fracture_99/NormalVec');

Frac(1).Radius = 20;
Frac(2).Radius = 20;
Frac(3).Radius = 5;
Frac(4).Radius = 5;

Frac(1).Center = [7.5; 0; 7.5];
Frac(2).Center = [7.5; 0; -7.5];
Frac(3).Center = [7.5; 0; 7.5];
Frac(4).Center = [7.5; 0; -7.5];

Frac(1).NormalVec = [-1.0; 0.; 0.5];
Frac(2).NormalVec = [1.0; 0.; 0.5];
Frac(3).NormalVec = [0.0; 1; 0.];
Frac(4).NormalVec = [0.0; 1; 0.];

for i = 1:4

    Frac(i).Conductivity = ((1.0e-3.^3.0) / 12.0);

    Frac(i).Verts3DTruncated = zeros(8, 3);
    Frac(i).NumVertsTruncated = 4;
    Frac(i).ConnectModelSurf = zeros(6, 1);

    Frac(i).NormalVec = Frac(i).NormalVec ./ norm(Frac(i).NormalVec);

    P_vec = [-Frac(i).NormalVec(2); Frac(i).NormalVec(1); 0];

    angle_degree = 45.0;

    Frac(i).Verts3D(1, :) = Quaternion_Rotation(angle_degree, Frac(i).NormalVec(1), Frac(i).NormalVec(2), Frac(i).NormalVec(3), P_vec(1), P_vec(2), P_vec(3));

    Frac(i).Verts3D(1, :) = Frac(i).Verts3D(1, :) ./ norm(Frac(i).Verts3D(1, :));

    Frac(i).Verts3D(1, :) = Frac(i).Verts3D(1, :) .* Frac(i).Radius;

    Frac(i).Verts3D(2, :) = cross(Frac(i).Verts3D(1, :), Frac(i).NormalVec);
    Frac(i).Verts3D(2, :) = Frac(i).Verts3D(2, :) ./ norm(Frac(i).Verts3D(2, :));
    Frac(i).Verts3D(2, :) = Frac(i).Verts3D(2, :) .* Frac(i).Radius;

    Frac(i).Verts3D(3, :) = -Frac(i).Verts3D(1, :) + Frac(i).Center';
    Frac(i).Verts3D(4, :) = -Frac(i).Verts3D(2, :) + Frac(i).Center';

    Frac(i).Verts3D(1, :) = Frac(i).Verts3D(1, :) + Frac(i).Center';
    Frac(i).Verts3D(2, :) = Frac(i).Verts3D(2, :) + Frac(i).Center';

end

figure(1)
view(3)
pbaspect([1, 1, 1]); hold on
xlabel('x')
ylabel('y')
zlabel('z'); hold on
L = 0.5 * 30;
cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; -L, L, -L; -L, L, L];
patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4; 5 6 7 8; 9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on
patch('Vertices', Frac(1).Verts3D, 'Faces', [1, 2, 3, 4], 'FaceVertexCData', Frac(1).Verts3D(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 1); hold on;
patch('Vertices', Frac(2).Verts3D, 'Faces', [1, 2, 3, 4], 'FaceVertexCData', Frac(2).Verts3D(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 1); hold on;
patch('Vertices', Frac(3).Verts3D, 'Faces', [1, 2, 3, 4], 'FaceVertexCData', Frac(3).Verts3D(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 1); hold on;
patch('Vertices', Frac(4).Verts3D, 'Faces', [1, 2, 3, 4], 'FaceVertexCData', Frac(4).Verts3D(:, 3), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 1); hold on;

ui = 1.1;
xlim([-ui * L, ui * L])
ylim([-ui * L, ui * L])
zlim([-ui * L, ui * L])

h5create([currentPath, '/Fractures.h5'], '/L', size(L));
h5write([currentPath, '/Fractures.h5'], '/L', L * 2)

h5create([currentPath, '/Fractures.h5'], '/NumFractures', size(L));
h5write([currentPath, '/Fractures.h5'], '/NumFractures', size(Frac, 2));

for i = 1:size(Frac, 2)
    h5create([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/Conductivity'], size(Frac(i).Conductivity));
    h5create([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/Verts3D'], size(Frac(i).Verts3D));
    h5create([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/Center'], size(Frac(i).Center));
    h5create([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/Verts3DTruncated'], size(Frac(i).Verts3DTruncated));
    h5create([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/NumVertsTruncated'], size(Frac(i).NumVertsTruncated));
    h5create([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/Radius'], size(Frac(i).Radius));
    h5create([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/ConnectModelSurf'], size(Frac(i).ConnectModelSurf));
    h5create([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/NormalVec'], size(Frac(i).NormalVec));

    h5write([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/Conductivity'], (Frac(i).Conductivity));
    h5write([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/Verts3D'], (Frac(i).Verts3D));
    h5write([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/Center'], (Frac(i).Center));
    h5write([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/Verts3DTruncated'], (Frac(i).Verts3DTruncated));
    h5write([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/NumVertsTruncated'], (Frac(i).NumVertsTruncated));
    h5write([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/Radius'], (Frac(i).Radius));
    h5write([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/ConnectModelSurf'], (Frac(i).ConnectModelSurf));
    h5write([currentPath, '/Fractures.h5'], ['/Fracture_', num2str(i), '/NormalVec'], (Frac(i).NormalVec));

end
