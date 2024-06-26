clc
clear all
close all

syms x
kappaVec = [0; 0; 0; 0; 0; 0; 15; 0; 2.5; 17.5; 5.0; 9.5; 13.5; 20; 10; 14; 8; 11];

SizePara = [3, 7.5, 0, 0, 0;
            3, 7.5, 0, 0, 0;
            3, 7.5, 0, 0, 0;
            3, 7.5, 0, 0, 0;
            3, 7.5, 0, 0, 0;
            3, 7.5, 0, 0, 0;
            0, 1.5, 1, 15, 0;
            2, 1, 15, 0, 0;
            1, 8, 5.5, 1, 15;
            0, 2.1, 1, 15, 0;
            0, 3.0, 1, 15, 0; 
            0, 1.3, 1, 100, 0;
            0, 1.7, 1, 30, 0;
            1, 8, 9.5, 1, 15;
            0, 1.6, 1, 25, 0;
            0, 1.8, 1, 35, 0;
            0, 2.3, 1, 55, 0;
            0, 1.9, 1, 200, 0];
        
Vex = zeros(size(kappaVec, 1), 1);
meanR = Vex;
L_charact = Vex;

for i = 1:size(kappaVec, 1)
    PreFactor = 0;
    kappa = kappaVec(i);

    if (kappa == 0)
        PreFactor = 0.5;
    else
        PreFactor = 2 / pi * pi / 4 * ...
            2 / (sinh(kappaVec(i)) ^ 2) * (besseli(0, 2 * kappaVec(i)) - ...
            1 / kappaVec(i) * besseli(1, 2 * kappaVec(i)));
    end

    Mode = SizePara(i, 1);

    if (Mode == 0)
        alpha_ = SizePara(i, 2);

        max_ = SizePara(i, 3);
        min_ = SizePara(i, 4);

        f = (1 - alpha_) * x ^ (-alpha_) / ...
            (max_ ^ (1 - alpha_) - min_ ^ (1 - alpha_));

        R_cubic = eval(int(x ^ 3 * f, x, min_, max_));

        Vex(i) = [PreFactor * 8 * 2 ^ 0.5 * R_cubic];

        meanR(i) = [eval(int(x ^ 1 * f, x, min_, max_))];

        meanR_cubic = eval(int(x^3 * f, x, min_, max_));
        meanR_6 = eval(int(x^6 * f, x, min_, max_)); 
        L_charact(i) = (2^3 * meanR_6 / (meanR_cubic * 2^1.5))^(1/3);

    elseif (Mode == 1)
        max_ = SizePara(i, 5);
        min_ = SizePara(i, 4);

        ex = SizePara(i, 2);
        dx = SizePara(i, 3);

        mu = log(ex * ex / ((dx + ex * ex) ^ 0.5));
        sigma = log(1 + (dx) / (ex * ex)) ^ 0.5;

        f = 1.0 / (x * (2 * pi) ^ 0.5 * sigma) * ...
            exp(-1.0 / (2 * sigma ^ 2) * (log(x) - mu) ^ 2);

        F_b = eval(int(f, x, -inf, max_));
        F_a = eval(int(f, x, -inf, min_));

        S = [15];

        S(1) = int(x * f, x, min_, max_) / (F_b - F_a);

        meanR(i) = [S];

        S(1) = int(x ^ 3 * f, x, min_, max_) / (F_b - F_a);

        R_cubic = S;

        Vex(i) = [PreFactor * 8 * 2 ^ 0.5 * R_cubic];

        meanR_cubic = eval(int(x^3 * f, x, min_, max_)) / (F_b - F_a);
        meanR_6 = eval(int(x^6 * f, x, min_, max_)) / (F_b - F_a); 

        L_charact(i) = (2^3 * meanR_6 / (meanR_cubic * 2^1.5))^(1/3);

    elseif (Mode == 2)
        f = 1 / (SizePara(i, 3) - SizePara(i, 2));

        R_cubic = eval(int(x ^ 3 * f, x, SizePara(i, 2), SizePara(i, 3)));

        Vex(i) = [PreFactor * 8 * 2 ^ 0.5 * R_cubic];

        meanR(i) = [eval(int(x ^ 1 * f, x, SizePara(i, 2), SizePara(i, 3)))];
        
        min_ = SizePara(i, 2);
        max_ = SizePara(i, 3);

        meanR_cubic = eval(int(x^3 * f, x, min_, max_));
        meanR_6 = eval(int(x^6 * f, x, min_, max_)); 
        L_charact(i) = (2^3 * meanR_6 / (meanR_cubic * 2^1.5))^(1/3);

    elseif (Mode == 3)
        Vex(i) = [PreFactor * 8 * 2 ^ 0.5 * SizePara(i, 2) .^ 3];

        meanR(i) = [SizePara(i, 2)];
        
        L_charact(i) = (SizePara(i, 2) ^ 2 * 2)^0.5;
    end

end
