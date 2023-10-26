clc
clear all
close all

rows = 1;
cols = 3;

currentPath = fileparts(mfilename('fullpath'));

pathYsd = '/home/tingchangyin/cuDFNsys/CompareCPUGPU/CPU_check_DFNgen_Nproc_10';
filename_rr = [pathYsd, '/DFN_gen_countTime.h5'];

pathYsk = '/home/tingchangyin/cuDFNsys/CompareCPUGPU/GPU_check_DFNgen';
filename_kk = [pathYsk, '/DFN_gen_countTime.h5'];

DensityIncreament = h5read(filename_rr, ['/DensityIncrement']);
DensityIncreament = double(DensityIncreament);
EndLoopStep = double(h5read(filename_rr, ['/EndLoopStep']));

CPU_n_I = [];
CPU_mean_DFN_gen_time = [];
CPU_sum_DFN_gen_time = [];

GPU_n_I = [];
GPU_mean_DFN_gen_time = [];
GPU_sum_DFN_gen_time = [];

for i = [1:1:EndLoopStep]
    DFNgenTime = h5read(filename_rr, ['/Step_', num2str(i, '%05d'), '/DFNgenTime']);
    n_I = h5read(filename_rr, ['/Step_', num2str(i, '%05d'), '/n_I']);
    DFNgenTimeall = h5read(filename_rr, ['/Step_', num2str(i, '%05d'), '/SumOfTime']);

    CPU_n_I = [CPU_n_I; 2 * mean(n_I)];
    CPU_mean_DFN_gen_time = [CPU_mean_DFN_gen_time; mean(DFNgenTime)];
    CPU_sum_DFN_gen_time = [CPU_sum_DFN_gen_time; DFNgenTimeall];

    DFNgenTime2 = h5read(filename_kk, ['/Step_', num2str(i, '%05d'), '/DFNgenTime']);
    n_I2 = h5read(filename_kk, ['/Step_', num2str(i, '%05d'), '/n_I']);
    DFNgenTimeall2 = h5read(filename_kk, ['/Step_', num2str(i, '%05d'), '/SumOfTime']);

    GPU_n_I = [GPU_n_I; 2 * mean(n_I2)];
    GPU_mean_DFN_gen_time = [GPU_mean_DFN_gen_time; mean(DFNgenTime2)];
    GPU_sum_DFN_gen_time = [GPU_sum_DFN_gen_time; DFNgenTimeall2];
end

figure(1)
subplot(rows, cols, 1)
title('(a)'); hold on
a1 = plot([1:1:EndLoopStep]' * DensityIncreament(1), CPU_n_I, 'r-+'); hold on
a2 = plot([1:1:EndLoopStep]' * DensityIncreament(1), GPU_n_I, 'b-o'); hold on
legend([a1 a2], 'CPU', 'GPU');
xlabel('$NF$', 'interpreter', 'latex')
ylabel('$n_I$', 'interpreter', 'latex')

figure(1)
subplot(rows, cols, 2)
title('(b)'); hold on
s1 = plot([1:1:EndLoopStep]' .* DensityIncreament, CPU_mean_DFN_gen_time, 'r-+'); hold on
s2 = plot([1:1:EndLoopStep]' .* DensityIncreament, GPU_mean_DFN_gen_time, 'b-o'); hold on
legend([s1 s2], 'CPU', 'GPU');
xlabel('$NF$', 'interpreter', 'latex')
ylabel('Second', 'interpreter', 'latex')

figure(1)
subplot(rows, cols, 3)
title('(c)'); hold on
s3 = plot([1:1:EndLoopStep]' .* DensityIncreament, CPU_sum_DFN_gen_time, 'r-+'); hold on
s4 = plot([1:1:EndLoopStep]' .* DensityIncreament, GPU_sum_DFN_gen_time, 'b-o'); hold on
legend([s3 s4], 'CPU', 'GPU');
xlabel('$NF$', 'interpreter', 'latex')
ylabel('Second', 'interpreter', 'latex')

for i = 1:3
    figure(1)
    subplot(rows, cols, i)
    hold on; set(gca,'FontSize',14); 
end



% ---------------------MHFEM triplets count time
rows = 1;
cols = 2;
pathYer = '/home/tingchangyin/cuDFNsys/CompareCPUGPU/CPU_check_MHFEM_Triplet_Nproc_10';
filename_qq = [pathYer, '/MHFEM_Triplet_countTime.h5'];

pathYwq = '/home/tingchangyin/cuDFNsys/CompareCPUGPU/GPU_check_MHFEM_Triplet';
filename_zz = [pathYwq, '/MHFEM_Triplet_countTime.h5'];

MC_times = h5read(filename_qq, '/MonteCarloTimes');

EndLoopStep = 20;

CPU_element_NUM = [];
CPU_Q_in = [];
CPU_Q_out = [];
CPU_TripletTime = [];

GPU_element_NUM = [];
GPU_Q_in = [];
GPU_Q_out = [];
GPU_TripletTime = [];

for i = 1:EndLoopStep

    CPU_element_NUM = [CPU_element_NUM; mean(h5read(filename_qq, ['/Step_', num2str(i, '%05d'), '/NUMeles']))];
    CPU_Q_in = [CPU_Q_in; mean(h5read(filename_qq, ['/Step_', num2str(i, '%05d'), '/Q_in']))];
    CPU_Q_out = [CPU_Q_out; mean(h5read(filename_qq, ['/Step_', num2str(i, '%05d'), '/Q_out']))];
    CPU_TripletTime = [CPU_TripletTime; mean(h5read(filename_qq, ['/Step_', num2str(i, '%05d'), '/TripletTime']))];

    GPU_element_NUM = [GPU_element_NUM; mean(h5read(filename_zz, ['/Step_', num2str(i, '%05d'), '/NUMeles']))];
    GPU_Q_in = [GPU_Q_in; mean(h5read(filename_zz, ['/Step_', num2str(i, '%05d'), '/Q_in']))];
    % GPU_Q_out = [GPU_Q_out; mean(h5read(filename_zz, ['/Step_', num2str(i, '%05d'), '/Q_out']))];
    GPU_TripletTime = [GPU_TripletTime; mean(h5read(filename_zz, ['/Step_', num2str(i, '%05d'), '/TripletTime']))];

end

[uniquevalue1, indexOfunique1] = unique(CPU_element_NUM, 'rows');

CPU_element_NUM = CPU_element_NUM(indexOfunique1, 1);
CPU_Q_in = CPU_Q_in(indexOfunique1, 1);
CPU_Q_out = CPU_Q_out(indexOfunique1, 1);
CPU_TripletTime = CPU_TripletTime(indexOfunique1, 1);

[uniquevalue2, indexOfunique2] = unique(GPU_element_NUM, 'rows');

GPU_element_NUM = GPU_element_NUM(indexOfunique2, 1);
GPU_Q_in = GPU_Q_in(indexOfunique2, 1);
% GPU_Q_out = GPU_Q_out(indexOfunique2, 1);
GPU_TripletTime = GPU_TripletTime(indexOfunique2, 1);

figure(2)
subplot(rows, cols, 1)
title('(a)'); hold on
w1 = plot(log(CPU_element_NUM), CPU_Q_in, 'r-+'); hold on
w2 = plot(log(GPU_element_NUM), GPU_Q_in, 'b-o'); hold on
legend([w1 w2], 'CPU', 'GPU');
xlabel('$log(NT)$', 'interpreter', 'latex')
ylabel('$Q_h \left[L^3T^{-1}\right]$', 'interpreter', 'latex')
ylim([0, mean(CPU_Q_in) * 2])

figure(2)
subplot(rows, cols, 2)
title('(b)'); hold on
k1 = plot(log(CPU_element_NUM), CPU_TripletTime, 'r-+'); hold on
k2 = plot(log(GPU_element_NUM), GPU_TripletTime, 'b-o'); hold on
legend([k1 k2], 'CPU', 'GPU');
xlabel('$log(NT)$', 'interpreter', 'latex')
ylabel('Second', 'interpreter', 'latex')

for i = 1:2
    figure(2)
    subplot(rows, cols, i)
    hold on; set(gca,'FontSize',14); 
end