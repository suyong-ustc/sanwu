function [x, u_mean, u_std] = PolyDisplacementByDIC(L, N, M)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    %% 读取DIC计算数据
    
    % 读取数据
    [x, ~, u, ~, zncc, ~] = ReadPolyDICData(L, N, M);
    
    invalid = zncc < 0.8;
    num_invalid = sum(sum(invalid));
    if num_invalid ~= 0
        disp("存在无效计算点！");
    end

    % 提取位移
    u_mean = mean(u);
    u_std = std(u);

    if L == 4
        x = x(1,:) - 40;
    else
        x = x(1,:) - 150;
    end

end