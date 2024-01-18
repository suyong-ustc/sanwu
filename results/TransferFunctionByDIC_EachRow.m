function fit_coeff = TransferFunctionByDIC_EachRow(T, M, U, N)
%AnalyzeTransferFunction 分析不满子区的传递函数
%   读取DIC计算结果，分析周期为T，子区大小M，空缺U，形函数为N情况下
%   正弦型变形的振幅衰减以及相位改变

    %% 读取DIC计算数据
    
    % 读取数据
    [x, ~, u, ~, zncc, ~] = ReadSinDICData(T, M, U, N);

    [r, ~] = size(x);
    
    if ~isempty(zncc(zncc<0.6))
        disp("存在无效计算点！");
    end
  
        
    %% 拟合正弦

    % 提取振幅
    t = x(1,:);
    
    % 待拟合函数
    fit = @(b,t) b(1) .* (sin(2*pi*t./T + b(2)));
   
    % 拟合系数
    fit_coeff = zeros(r, 2);

    % 拟合设置
    options = optimoptions('lsqcurvefit', 'Display','off');
    lb = [-1, -pi];
    ub = [2, pi];
    
    for j = 1 : r
        % 第 j 行位移
        ur = u(j,:);
    
        % 拟合振幅与相位
        s = lsqcurvefit(fit, [1; 0], t', ur', lb, ub, options);
    
        % 存储拟合结果
        fit_coeff(j,:) = s;
    end

end