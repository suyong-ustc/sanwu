function [H0, H1, H2] = TransferFunctionByTheory(T, M, U)
%   根据理论估计周期为T，子区大小M，空缺U情况下
%   正弦型变形的振幅衰减以及相位改变

    % 角频率
    omega = 2 * pi / T;
    
    % 无量纲化角频率
    mu = 0.5 * M * omega;
    
    % 无量纲化空缺度
    alpha = U / M;
    
    % 参数
    a = mu .* (1 - alpha);
    b = exp(-1j * mu .* alpha);
    
    s = sinc(a/pi);
    c = cos(a);
    
    % 零阶形函数传递函数
    H0 = b .* s;

    % 一阶形函数传递函数
    T1 = 1j * 3.0 * alpha ./ (a .* (1-alpha));
    H1 = H0 + b .* T1 .* (s - c);

    % 二阶形函数传递函数
    T2 = 2.5 * (3 * (alpha./(1-alpha)).^2 - 1);
    H2 = H1 + b .* T2 .* (s + 3 * (c - s) .* a.^-2);

end