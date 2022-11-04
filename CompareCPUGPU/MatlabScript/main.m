clc
clear all
close all

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
subplot(1, 3, 1)
title('Number of intersections per fracture'); hold on
a1 = plot([1:1:EndLoopStep]' * DensityIncreament(1), CPU_n_I, 'r-+'); hold on
a2 = plot([1:1:EndLoopStep]' * DensityIncreament(1), GPU_n_I, 'b-o'); hold on
legend([a1 a2], 'CPU', 'GPU');
xlabel('NF')
ylabel('n_I')

figure(1)
subplot(1, 3, 2)
title('Mean running time (s)'); hold on
s1 = plot([1:1:EndLoopStep]' .* DensityIncreament, CPU_mean_DFN_gen_time, 'r-+'); hold on
s2 = plot([1:1:EndLoopStep]' .* DensityIncreament, GPU_mean_DFN_gen_time, 'b-o'); hold on
legend([s1 s2], 'CPU', 'GPU');
xlabel('NF')
ylabel('t (s)')

figure(1)
subplot(1, 3, 3)
title('Sum of running time (s)'); hold on
s3 = plot([1:1:EndLoopStep]' .* DensityIncreament, CPU_sum_DFN_gen_time, 'r-+'); hold on
s4 = plot([1:1:EndLoopStep]' .* DensityIncreament, GPU_sum_DFN_gen_time, 'b-o'); hold on
legend([s3 s4], 'CPU', 'GPU');
xlabel('NF')
ylabel('t (s)')
