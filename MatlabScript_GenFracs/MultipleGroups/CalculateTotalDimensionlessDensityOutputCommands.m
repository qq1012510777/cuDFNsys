clc
clear all
close all

L = 250;
rhoc = 2.31;
phi_0 = [0, 0] * pi / 180;
theta_0 = [90, 0] * pi / 180;

kappa = [22, 14];
Weight = [0.3, 0.7];
size_para = [0, 2.4, 1, 25, 0;
    0, 1.3, 1, 50, 0];
Gamma = [1.2e-8, 3.3e-7];
Beta = [0.24, 0.31];

NumFractureSets = length(kappa);

Vex_para1 = 0; Vex_para2 = 0;
Vex = 0;
for i = 1:NumFractureSets
    theta_0_1 = theta_0(i);
    phi_0_1 = phi_0(i);
    kappa1 = kappa(i);
    W_1 = Weight(i);

    for j = 1:NumFractureSets
        theta_0_2 = theta_0(j);
        phi_0_2 = phi_0(j);
        kappa2 = kappa(j);
        W_2 = Weight(j);

        % Fisher
        fisher_axial = @(theta, phi, theta0, phi0, kappa) ...
            (1/(2*pi)) .* (kappa .* sin(theta) ./ sinh(kappa)) .* ...
            cosh(kappa .* (cos(theta0).*cos(theta) + ...
                   sin(theta0).*sin(theta).*cos(phi - phi0)));
        sinGamma = @(theta1, phi1, theta2, phi2) ...
            sqrt(1 - (cos(theta1).*cos(theta2) + ...
              sin(theta1).*sin(theta2).*cos(phi1 - phi2)).^2);

        f4 = @(theta1, phi1, theta2, phi2) ...
            sinGamma(theta1, phi1, theta2, phi2) .* ...
            fisher_axial(theta1, phi1, theta_0_1, phi_0_1, kappa1) .* ...
            fisher_axial(theta2, phi2, theta_0_2, phi_0_2, kappa2);
        q4 = integralN(f4,0, 0.5*pi,0,pi*2,0,0.5*pi,0,2*pi,'AbsTol',1e-5,'RelTol',1e-3);
        q4 = real(q4);
        Vex_para1 = Vex_para1 + W_1 * W_2 * q4;
    end

    
    if size_para(i, 1) == 0
        Vex_para2 = Vex_para2 + W_1 * 8*sqrt(2) * truncated_powerlaw_moment(3, size_para(i,2), size_para(i,3), size_para(i,4));
    end
end
Vex = 2/pi * Vex_para1 * Vex_para2;

NF_c = rhoc * L^3 / Vex;

NF_range = linspace(0, 4*NF_c, 31);

NF_range = round(NF_range(2:end));
IncreNF = NF_range(2) - NF_range(1);

disp([num2str(L), ' 1 1 1 ', num2str(NumFractureSets), ' \'])
for i = 1:NumFractureSets
    theta_0_1 = theta_0(i);
    phi_0_1 = phi_0(i);
    W_1 = Weight(i);
    disp([num2str(kappa(i)), ' ', num2str(sph2cart_unit(theta_0_1, phi_0_1)'), ' ', num2str(round(NF_range(1) * W_1)), ' ', ...
        num2str(round(IncreNF * W_1)), ' ', num2str(length(NF_range)-1), ' ', ...
        num2str(size_para(i,1)), ' ', num2str(size_para(i, 2:end)), ...
        ' ', num2str(Beta(i)), ' ', num2str(Gamma(i)), ' \'])
end

%--------------functions here
function Mk = truncated_powerlaw_moment(k, alpha, Rm, RM)
    % truncated_powerlaw_moment - kth moment of a truncated power-law distribution
    %
    % Mk    : kth moment E[R^k]
    % k     : moment order (can be scalar or array)
    % alpha : power-law exponent
    % Rm    : minimum R
    % RM    : maximum R

    if abs(alpha - 1) < 1e-12
        error('alpha = 1 is a special case (logarithmic), handle separately.');
    end
    
    % Normalization constant
    C = (1 - alpha) / (RM^(1 - alpha) - Rm^(1 - alpha));
    
    % Moment formula (valid if k - alpha + 1 ≠ 0)
    Mk = C .* ( (RM.^(k - alpha + 1) - Rm.^(k - alpha + 1)) ./ (k - alpha + 1) );
end
function v = sph2cart_unit(theta0, phi0)
    % Convert spherical angles (theta0, phi0) to a normalized 3D vector
    
    x = sin(theta0) * cos(phi0);
    y = sin(theta0) * sin(phi0);
    z = cos(theta0);
    
    v = [x; y; z];  % 3x1 vector, unit length
end