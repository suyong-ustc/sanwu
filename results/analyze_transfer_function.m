clc;
clear;
close all;

%% 参数
T = 50;     % 正弦周期
M = 29;     % 子区尺寸
N = 1;      % 形函数阶数

ShowSinDICData(T,M,18,1);

U = 0:1:20;
a_mean = zeros(length(U),1);
a_std = zeros(length(U),1);
phi_mean = zeros(length(U),1);
phi_std = zeros(length(U),1);
for i = 1 : length(U)
    [a_mean(i), phi_mean(i)] = TransferFunctionByDIC_Mean(T, M, U(i), N);
    % [a_mean(i), a_std(i), phi_mean(i), phi_std(i)] = TransferFunctionByDIC_EachRow(T, M, U(i), N);
end


%% 理论估计
% 振幅的估计

[H0, H1, H2] = TransferFunctionByTheory(T, M, U);


%% 理论估计
% 振幅的估计

figure;

subplot(1, 2, 1);
hold on;
plot(U, a_mean, '*');
if N == 0
    plot(U, abs(H0), 'r');
elseif N == 1
    plot(U, abs(H1), 'r');
elseif N == 2
    plot(U, abs(H2), 'r');
end
hold off;

subplot(1, 2, 2);
hold on;
plot(U, phi_mean, '*');
if N == 0
    plot(U, angle(H0), 'r');
elseif N == 1
    plot(U, angle(H1), 'r');
elseif N == 2
    plot(U, angle(H2), 'r');
end
hold off;