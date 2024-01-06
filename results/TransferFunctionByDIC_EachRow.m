function [a_mean, a_std, phi_mean, phi_std] = TransferFunctionByDIC_EachRow(T, M, U, N)
%AnalyzeTransferFunction 分析不满子区的传递函数
%   读取DIC计算结果，分析周期为T，子区大小M，空缺U，形函数为N情况下
%   正弦型变形的振幅衰减以及相位改变

    %% 读取DIC计算数据
    
    % 读取数据
    [x, ~, u, ~, ~, s] = ReadSinDICData(T, M, U, N);

    [r, ~] = size(x);
    
    if ~isempty(s(s~=1))
        disp("存在无效计算点！");
    end
    
    u(s~=1) = nan;
  
        
    %% 拟合正弦

    % 提取振幅
    t = x(1,:);
    
    % 待拟合函数
    fit = @(b,t) b(1) .* (sin(2*pi*t./T + b(2)));
   
    % 拟合系数
    fit_coeff = zeros(r, 2);
    
    for j = 1 : r
        % 第 j 行位移
        ur = u(j,:);

        if sum(isnan(ur)) ~= 0
            fit_coeff(j,:) = NaN;
            continue;
        end
 
        % 位移场正弦拟合
        %fcn = @(b) sum((fit(b,t) - ur).^2);
    
        % 拟合振幅与相位
        s = lsqcurvefit(fit, [1; 0], t', ur');
    
        % 存储拟合结果
        fit_coeff(j,:) = s;
    end
    
    fit_amplitude = fit_coeff(:,1);
    a_mean = mean(fit_amplitude,"omitnan");
    a_std = std(fit_amplitude,"omitnan");
    
    fit_phase = fit_coeff(:,2);
    phi_mean = mean(fit_phase,"omitnan");
    phi_std = std(fit_phase,"omitnan");

end