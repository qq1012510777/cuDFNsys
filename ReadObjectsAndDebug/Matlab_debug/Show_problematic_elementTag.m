clc;
close all;
clear all;
load('MHFEM_1.mat');
L = 0.5 * 100;
cube_frame = [-L, -L, L; -L, L, L; L, L, L; L -L, L; -L, -L, -L; -L, L, -L; L, L, -L; L -L, -L; -L, L, L; -L, L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; L, -L, -L; L, -L, L; L, -L, L; L, -L, -L; -L, -L, -L; -L, -L, L; L, L, L; L, L, -L; -L, L, -L; -L, L, L];
figure(1); view(3); title('DFN flow (mhfem)'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on
patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4; 5 6 7 8; 9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on

xlim([-1.1 * L, 1.1 * L]);
ylim([-1.1 * L, 1.1 * L]);
zlim([-1.1 * L, 1.1 * L]);
hold on
patch('Vertices', coordinate_3D, 'Faces', element_3D, 'FaceVertexCData', pressure_eles, 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 1); colorbar; view(3); hold on
caxis([20, 100]);

CenterELE = zeros(size(element_3D, 1), 3);
CenterELE(:, 1) = 1/3 * (coordinate_3D(element_3D(:, 1), 1) + coordinate_3D(element_3D(:, 2), 1) + coordinate_3D(element_3D(:, 3), 1));
CenterELE(:, 2) = 1/3 * (coordinate_3D(element_3D(:, 1), 2) + coordinate_3D(element_3D(:, 2), 2) + coordinate_3D(element_3D(:, 3), 2));
CenterELE(:, 3) = 1/3 * (coordinate_3D(element_3D(:, 1), 3) + coordinate_3D(element_3D(:, 2), 3) + coordinate_3D(element_3D(:, 3), 3));
quiver3(CenterELE(:, 1), CenterELE(:, 2), CenterELE(:, 3), velocity_center_grid(:, 1), velocity_center_grid(:, 2), velocity_center_grid(:, 3), 4, 'LineWidth', 1.5, 'color', 'r');

close 1

S = load('DFN_mesh_1.mat');
NumFracs = S.element_Frac_Tag(end, 1);

Frac = [];

for i = 1:NumFracs
    F = h5read('Fractures_II.h5', ['/Fracture_', num2str(i), '/Verts3D']);
    Frac(:, :, i) = F;
    clear F
end

dist = [];

for i = 1:size(element_3D, 1)
    FracID = S.element_Frac_Tag(i, 1);
    Plane = Frac([1, 2, 3], :, FracID);

    for j = 1:3
        nodeID = element_3D(i, j);

        pnt = coordinate_3D(nodeID, :);

        x1 = Plane(0 + 1, 1);
        x2 = Plane(1 + 1, 1);
        x3 = Plane(2 + 1, 1);
        y1 = Plane(0 + 1, 2);
        y2 = Plane(1 + 1, 2);
        y3 = Plane(2 + 1, 2);
        z1 = Plane(0 + 1, 3);
        z2 = Plane(1 + 1, 3);
        z3 = Plane(2 + 1, 3);
        a1 = x2 - x1;
        b1 = y2 - y1;
        c1 = z2 - z1;
        a2 = x3 - x1;
        b2 = y3 - y1;
        c2 = z3 - z1;

        a = b1 * c2 - b2 * c1;
        b = a2 * c1 - a1 * c2;
        c = a1 * b2 - b1 * a2;
        d = (-a * x1 - b * y1 - c * z1);

        d = abs((a * pnt(1) + b * pnt(2) + c * pnt(3) + d));

        e = (a * a + b * b + c * c) ^ 0.5;

        dist(i, j) = d / e;
    end

end

% a = find(S.element_Frac_Tag == 12);
% figure(2); view(3);
% patch('Vertices', coordinate_3D, 'Faces', element_3D(a, :), 'FaceVertexCData', pressure_eles(a, :), 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 0); colorbar; view(3); hold on
% 
% c = 1828;
% patch('Vertices', coordinate_3D, 'Faces', element_3D(c, :), 'FaceVertexCData', pressure_eles(c, :), 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 1, 'edgecolor', 'r', 'linewidth', 2.5); colorbar; view(3); hold on
% 
% scatter3(coordinate_3D(351, 1), coordinate_3D(351, 2), coordinate_3D(351, 3), 'o', 'k');
% 
% figure(3)
% patch('Vertices', coordinate_3D, 'Faces', element_3D(c, :), 'FaceVertexCData', pressure_eles(c, :), 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 0, 'edgecolor', 'k', 'linewidth', 1.0); colorbar; view(3); hold on
% 
% coordinate = [-2.98168, 10.2286, -41.6272
% -2.90953, 10.1941, -41.6421
% -2.97118, 10.2421, -41.6369
% 
% ];
% hold on
% patch('Vertices', coordinate, 'Faces', [1 2 3], 'FaceVertexCData', [1], 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 0, 'edgecolor', 'r', 'linewidth', 1.0); colorbar; view(3); hold on
