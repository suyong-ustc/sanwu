function [a, phi] = TransferFunctionByDIC_Mean(T, M, U, N)
%AnalyzeTransferFunction 分析不满子区的传递函数
%   读取DIC计算结果，分析周期为T，子区大小M，空缺U，形函数为N情况下
%   正弦型变形的振幅衰减以及相位改变

    disp(['T=',num2str(T),'; M=',num2str(M),'; U=',num2str(U),'; N=',num2str(N)]);

    %% 读取DIC计算数据
    
    % 读取数据
    [x, ~, u, ~, ~, s] = ReadSinDICData(T, M, U, N);
    
    if ~isempty(s(s~=1))
        disp("存在无效计算点！");
    end

    u(s~=1) = nan;
     
    %% 拟合正弦
    
    % 提取振幅
    t = x(1,:);
    
    % 待拟合函数
    fit = @(b,t) b(1) .* (sin(2*pi*t./T + b(2)));

    % 拟合振幅与相位
    [s,resnorm,residual,exitflag,output] = lsqcurvefit(fit, [1; 0], t, mean(u,"omitnan"));

    a = s(1);
    phi = s(2);

end