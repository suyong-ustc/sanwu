clc;
clear;
close all;

%% 参数
T = 60;     % 正弦周期
M = 29;     % 子区尺寸
N = 2;      % 形函数阶数

% ShowSinDICData(T,M,9,N);

%% 分析DIC计算结果
U = 0 : 1 : (M+1)/2;

transfer_func_by_row = zeros(4, length(U));
transfer_func_by_mean = zeros(2, length(U));


for i = 1 : length(U)

    disp(['T=',num2str(T),'; M=',num2str(M),'; U=',num2str(U(i)),'; N=',num2str(N)]);

    % 对DIC每列的平均位移进行拟合
    [transfer_func_by_mean(1,i), transfer_func_by_mean(2,i)] = TransferFunctionByDIC_Mean(T, M, U(i), N);

    % 每行进行拟合
    fit_coeff = TransferFunctionByDIC_EachRow(T, M, U(i), N);

    transfer_func_by_row(1,i) = mean(fit_coeff(:,1));
    transfer_func_by_row(2,i) = std(fit_coeff(:,1));
    transfer_func_by_row(3,i) = mean(fit_coeff(:,2));
    transfer_func_by_row(4,i) = std(fit_coeff(:,2));

end


%% 理论估计
% 振幅的估计

UU = 0 : 0.1 : 1+M/2;
[H0, H1, H2] = TransferFunctionByTheory(T, M, UU);


%% 绘图

figure;

subplot(1, 2, 1);
hold on;
errorbar(U, transfer_func_by_row(1,:), transfer_func_by_row(2,:));
plot(U, transfer_func_by_mean(1,:), 'k*');
if N == 0
    plot(UU, abs(H0), 'r');
elseif N == 1
    plot(UU, abs(H1), 'r');
elseif N == 2
    plot(UU, abs(H2), 'r');
end
hold off;

subplot(1, 2, 2);
hold on;
errorbar(U, transfer_func_by_row(3,:), transfer_func_by_row(4,:));
plot(U, transfer_func_by_mean(2,:), 'k*');
if N == 0
    plot(UU, angle(H0), 'r');
elseif N == 1
    plot(UU, angle(H1), 'r');
elseif N == 2
    plot(UU, angle(H2), 'r');
end
hold off;


%% 数据存储
data1 = [U', transfer_func_by_row'];
if N == 0
    data2 = [UU', abs(H0)', angle(H0)'];
elseif N == 1
    data2 = [UU', abs(H1)', angle(H1)'];
elseif N == 2
    data2 = [UU', abs(H2)', angle(H2)'];
end