clc
clear all
close all

syms x
% ----------sin gamma
kappaVec = [0; 15; 7.5; 2.5; 17.5; 5.0];
PreFactor = []; % sin gamma

for i = 1:size(kappaVec, 1)

    if (kappaVec(i) == 0)
        PreFactor(i, 1) = 1/2;
    else
        PreFactor(i, 1) = 2 / pi * pi / 4 * ...
            2 / (sinh(kappaVec(i)) ^ 2) * (besseli(0, 2 * kappaVec(i)) - ...
            1 / kappaVec(i) * besseli(1, 2 * kappaVec(i)));
    end

end

%---------- power laws
alphaVec = [1.5, 2.1, 3.0];
minRVec = [1, 1, 1];
maxRVec = [15, 15, 15];

Vex = [];
meanR = [];
m = 1;

for i = 1:size(kappaVec, 1)

    for j = 1:size(alphaVec, 2)

        alpha_ = alphaVec(j);

        max_ = maxRVec(j);
        min_ = minRVec(j);

        f = (1 - alpha_) * x ^ (-alpha_) / ...
            (max_ ^ (1 - alpha_) - min_ ^ (1 - alpha_));

        R_cubic = eval(int(x ^ 3 * f, x, min_, max_));

        Vex = [Vex; PreFactor(i) * 8 * 2 ^ 0.5 * R_cubic];
        
        meanR = [meanR; eval(int(x ^ 1 * f, x, min_, max_))];

        m = m + 1;
        clear f alpha_ max_ min_
    end

end

% --------monosize
%---------- power laws

RVec = [7.5];

Vex_mono = [];
meanR_mono = [];
m = 1;

for i = 1:size(RVec, 1)

    for j = 1:size(alphaVec, 1)

        R = RVec(j);
       
        Vex_mono = [Vex_mono; PreFactor(i) * 8 * 2 ^ 0.5 * R .^ 3];
        
        meanR_mono = [meanR_mono; R];

        m = m + 1;
    end

end

% ---------uniform
Rmin=[1];
Rmax=[15];
m = 1;

V_ex_ufm = [];
meanR_ufm = [];

for i = 1:size(kappaVec, 1)

    for j = 1:size(Rmin, 1)

        f = 1 / (Rmax(j) - Rmin(j));
        
        R_cubic = eval(int(x ^ 3 * f, x, Rmin(j), Rmax(j)));
        
        V_ex_ufm = [V_ex_ufm; PreFactor(i) * 8 * 2 ^ 0.5 * R_cubic];
        
        meanR_ufm = [meanR_ufm; eval(int(x ^ 1 * f, x, Rmin(j), Rmax(j)))];

        m = m + 1;
    end

end

% ---------lognormal
Rmin=[1];
Rmax=[15];
m_l = [8];
v_l = [5.5];

V_ex_lognormal = [];
meanR_lognormal = [];

m = 1;
for i = 1:size(kappaVec, 1)

    for j = 1:size(Rmin, 1)
        
        max_ = Rmin(j);
        min_ = Rmax(j);

        ex = m_l(j);
        dx = v_l(j);

        mu = log(ex * ex / ((dx + ex * ex)^0.5));
        sigma = log(1 + (dx) / (ex * ex))^0.5;
            
        f = 1.0 / (x * (2 * pi)^0.5 * sigma) * ...
                exp(-1.0 / (2 * sigma^2) * (log(x) - mu)^2);
            
        F_b = eval(int(f, x, -inf, max_));
        F_a = eval(int(f, x, -inf, min_));

        S = [15];

        S(1) = int(x * f, x, min_, max_) / (F_b - F_a);

        meanR_lognormal = [meanR_lognormal, ...
                            S];

        S(1) = int(x^3 * f, x, min_, max_) / (F_b - F_a);

        R_cubic = S;
                
        V_ex_lognormal = [V_ex_lognormal; PreFactor(i) * 8 * 2 ^ 0.5 * R_cubic];
        
        % meanR_lognormal = [meanR_lognormal; eval(int(x ^ 1 * f, x, min_, max_))];

        m = m + 1;
    end

end