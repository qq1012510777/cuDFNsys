clc
clear all
close all
syms x

fisher_no = 1;

kappa = [0.0, 2.5, 7.5];

Fai = [1, 0, 0, 0];
Sin_gamma_mean = [pi / 4, 0, 0, 0];
prefactor_Vex = [1/2, 0, 0, 0];

Psi = [0, 0, 0];

for i = 2:3
    kappa_1 = kappa(i);
    Fai(i) = 2 / (sinh(kappa_1)^2) * (besseli(0, 2 * kappa_1) - 1 / kappa_1 * besseli(1, 2 * kappa_1));
    Sin_gamma_mean(i) = Fai(i) * pi / 4;
    prefactor_Vex(i) = Sin_gamma_mean(i) * 2 / pi;

    Psi(i) = 3 / kappa_1^2 * (kappa_1 * coth(kappa_1) - 1);
end

sigma_exponent = [0.1, 0.5];

%--------power law
min_vec = 1;
max_vec = 100;

alpha_vec = [1.5, 2.0, 2.5];

mean_powerlaw = [];
E_R_3_powerlaw = [];
meansize_permeab_powerlaw = [];
E_R_6_powerlaw = [];
V_ex_powerlaw = [];
Lc_powerlaw = [];
Sigma_powerlaw = [];

expo_k = 4;
expo_j = 2;
expo_l = 3/4;
%-------------modelno
%-------------modelno
%-------------modelno
%-------------modelno
model_NO = [1:5:10];

for i = 1:size(alpha_vec, 2)

    for fisher_no = 1:1

        for k = 1:size(sigma_exponent, 2)

            alpha_ = alpha_vec(i);

            max_ = max_vec;
            min_ = min_vec;

            f = (1 - alpha_) * x^(-alpha_) / (max_^(1 - alpha_) - min_^(1 - alpha_));

            mean_powerlaw = [mean_powerlaw, eval(int(x * f, x, min_, max_))];

            E_R_3_powerlaw = [E_R_3_powerlaw, eval(int(x^3 * f, x, min_, max_))];

            V_ex_powerlaw = [V_ex_powerlaw, 8 * 2^0.5 * E_R_3_powerlaw(end) * prefactor_Vex(fisher_no)];

            meansize_permeab_powerlaw = [meansize_permeab_powerlaw, ...
                                        eval(int((2^0.5 * x)^expo_k * f, x, min_, max_))];

            % 2.31 / V_ex_powerlaw(i) * (L ^ 3) / 20
            E_R_6_powerlaw = [E_R_6_powerlaw, eval(int(x^6 * f, x, min_, max_))];

            % Lc_powerlaw(i) = ( (8 * (E_R_6_powerlaw(i) - E_R_3_powerlaw(i)^2)) / (2 * 2 ^ 0.5 * E_R_3_powerlaw(i))) ^(1/3);
            Lc_powerlaw = [Lc_powerlaw, (2^3 * E_R_6_powerlaw(end) / (E_R_3_powerlaw(end) * 2^1.5))^(1/3)];
            %Lc_powerlaw = [Lc_powerlaw, (2^1.5* (E_R_3_powerlaw(end)))^(1/3)];

            % fracture conductivity distribution
            Sigma_powerlaw = [Sigma_powerlaw, eval(int((1.0e-11) * (2^0.5 * x)^(expo_j) ...
                        * (x^(3 * sigma_exponent(k)))^(1) * f, min_, max_))];
        end

    end

end

mean_powerlaw = mean_powerlaw';
E_R_3_powerlaw = E_R_3_powerlaw';
meansize_permeab_powerlaw = meansize_permeab_powerlaw';
E_R_6_powerlaw = E_R_6_powerlaw';
V_ex_powerlaw = V_ex_powerlaw';
Lc_powerlaw = Lc_powerlaw';
Sigma_powerlaw = Sigma_powerlaw';
