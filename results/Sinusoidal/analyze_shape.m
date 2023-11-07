clc;
clear;
close all;

%% 参数
T = 40;     % 正弦周期
M = 19;     % 子区尺寸
D = 3;      % 空
N = 1;      % 形函数阶数

alpha = D / M;
omega = 2 * pi / T;
mu = 0.5 * M * omega;


%% 读取DIC计算结果
% 读取数据
prefix = ['T', num2str(T), 'N', num2str(N), 'D', num2str(D), 'M', num2str(M)];
dic_x = readmatrix([prefix, '__x.csv']);
dic_u = readmatrix([prefix, '__u.csv']);
dic_zncc = readmatrix([prefix, '__zncc.csv']);

invalid = (dic_zncc < 0.8);
disp(sum(sum(invalid)));

% 提取振幅
x = dic_x(1,:);
dic_u_mean = mean(dic_u);
dic_u_std = std(dic_u);


%% 拟合计算结果
% 第一行：amplitude average
% 第二行：amplitude standard deviation
% 第三行：phase average
% 第四行：phase standard deviation
dic_results = zeros(4,1);

fit = @(b,x) b(1) .* (sin(2*pi*x./b(2) + b(3)));

[r, c] = size(dic_u);

% 拟合
fit_results = zeros(3, c);
for j = 1 : c
    % 特定散斑场的计算位移场
    u = dic_u(:,j);

    % 位移场正弦拟合
    fcn = @(b) sum((fit(b,x) - u').^2);
    s = fminsearch(fcn, [1; T; 0]);

    % 存储拟合结果
    fit_results(:,j) = s;
end

% 均值与标准差
dic_results(1) = mean(fit_results(1,:));
dic_results(2) = std(fit_results(1,:));
dic_results(3) = mean(fit_results(3,:));
dic_results(4) = std(fit_results(3,:));



%% 理论估计
% 振幅的估计
s = 0;
c = 0;
for k = 0 : N
    fun_sin = @(z) ...
        sin(mu*(1-alpha)*z) .* legendreP(k,z);
    s = s + (k+0.5) * legendreP(k,alpha/(1-alpha)) * integral(fun_sin,-1,1);

    fun_cos = @(z) ...
        cos(mu*(1-alpha)*z) .* legendreP(k,z);
    c = c + (k+0.5) * legendreP(k,alpha/(1-alpha)) * integral(fun_cos,-1,1);
end

theory_amplitude = sqrt(s*s+c*c);

% 相位的估计
theory_phase = atan2(s,c)-mu*alpha;


figure;
hold on;
errorbar(x,dic_u_mean,dic_u_std);
fplot(@(x) theory_amplitude*sin(2*pi*x/T+theory_phase), [min(x) max(x)]);

%% 理论估计结果
% Xu的理论
% xu_ratio = xu_theoretical_amplitude(T,M,N);
% xu_u = xu_ratio * real_u;
% 
% 我的理论
% su_ratio = su_theoretical_amplitude(T,M,N);
% su_u = su_ratio * real_u;



% plot(x,xu_u,'x');
% plot(x,su_u,'o');

%data = [x',dic_u_mean',dic_u_std',xu_u',su_u'];
