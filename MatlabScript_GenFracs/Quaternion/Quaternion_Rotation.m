function f = Quaternion_Rotation(Angle_degree, axis_x, axis_y, axis_z, x, y, z)

    currentPath3 = fileparts(mfilename('fullpath'));
    addpath(genpath([currentPath3, '/Quaternion_cite']));

    Q = qGetRotQuaternion(Angle_degree * pi / 180, [axis_x, axis_y, axis_z]);

    f = qRotatePoint([x y z], Q);
    f = f';
end
