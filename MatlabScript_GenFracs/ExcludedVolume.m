clc
clear all
close all

syms x
% ----------sin gamma
kappaVec = [0];
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
alphaVec = [1.5];
minRVec = [1];
maxRVec = [15];

Vex = [];
meanR = [];
m = 1;

for i = 1:size(kappaVec, 1)

    for j = 1:size(alphaVec, 1)

        alpha_ = alphaVec(j);

        max_ = maxRVec(j);
        min_ = minRVec(j);

        f = (1 - alpha_) * x ^ (-alpha_) / ...
            (max_ ^ (1 - alpha_) - min_ ^ (1 - alpha_));

        R_cubic = eval(int(x ^ 3 * f, x, min_, max_));

        Vex = [Vex; PreFactor(i) * 8 * 2 ^ 0.5 * R_cubic];
        
        meanR = [meanR; eval(int(x ^ 1 * f, x, min_, max_))];

        m = m + 1;
    end

end
