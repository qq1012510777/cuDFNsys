clc;
close all;
clear all;
Trajectory = [	2.935369003019592692993455784745, -2.345783909633889052770427952055
	2.908520000000000216289208765374, -2.317600000000000104449782156735];

Nodes = [% ParticleID: 587, neighboring elementID: 10826, fracID: 58:
	1.687137718547782494482589754625, -3.103163235365630878703768757987
	2.956274559019663961123569606571, -2.319792750866472186288547163713
	1.868715640541841471744533009769, -2.090865832087615228118693266879
% ParticleID: 587, neighboring elementID: 10869, fracID: 58:
	1.687137718547782494482589754625, -3.103163235365630878703768757987
	1.868715640541841471744533009769, -2.090865832087615228118693266879
	0.529799318577021005616245474812, -2.601294285837996600463384311297
% ParticleID: 587, neighboring elementID: 10914, fracID: 58:
	0.418000878075931614485938325743, -3.886533719864846858627061010338
	1.687137718547782494482589754625, -3.103163235365630878703768757987
	0.529799318577021005616245474812, -2.601294285837996600463384311297
% ParticleID: 587, neighboring elementID: 10921, fracID: 58:
	1.868715640541841471744533009769, -2.090865832087615228118693266879
	2.956274559019663961123569606571, -2.319792750866472186288547163713
	2.627802336598911825404911724036, -1.589699784259942738628978986526
% ParticleID: 587, neighboring elementID: 10925, fracID: 58:
	2.956274559019663961123569606571, -2.319792750866472186288547163713
	3.487043747708560736953131709015, -1.659905673382432889084725502471
	2.627802336598911825404911724036, -1.589699784259942738628978986526
% ParticleID: 587, neighboring elementID: 10988, fracID: 58:
	1.687137718547782494482589754625, -3.103163235365630878703768757987
	0.418000878075931614485938325743, -3.886533719864846858627061010338
	1.325641191014091724298396002268, -4.347103134724975781466582702706
% ParticleID: 587, neighboring elementID: 10989, fracID: 58:
	2.140957875016876066354143404169, -3.333447942795725982279009258491
	1.687137718547782494482589754625, -3.103163235365630878703768757987
	1.325641191014091724298396002268, -4.347103134724975781466582702706
    ];

Triangles = [1:size(Nodes, 1)];
Triangles = reshape(Triangles, 3, size(Nodes, 1) / 3);
Triangles = Triangles';

Nodes_p = [1.687137718547782494482589754625, -3.103163235365630878703768757987
	2.140957875016876066354143404169, -3.333447942795725982279009258491
	2.956274559019663961123569606571, -2.319792750866472186288547163713];
Triangles_p = [1:size(Nodes_p, 1)];
Triangles_p = reshape(Triangles_p, 3, size(Nodes_p, 1) / 3);
Triangles_p = Triangles_p';

figure(1)
patch('Vertices', Nodes, 'Faces', Triangles, 'FaceVertexCData', zeros(size(Nodes, 1), 1), 'FaceColor', ...
    'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on
plot(Trajectory(:, 1), Trajectory(:, 2), '*-r'); hold on
patch('Vertices', Nodes_p, 'Faces', Triangles_p, 'FaceVertexCData', zeros(size(Nodes_p, 1), 1), 'FaceColor', ...
    'interp', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'b'); hold on