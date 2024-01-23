function [x, y, u, v, zncc, s] = ReadPolyDICData(L, N, M)
%ReadBoundaryDICData 读取多项式变形边缘的DIC分析结果
%   L : 多项式阶数
%   N : 形函数阶数
%   M : 子区尺寸
%   X : DIC计算结果的X坐标
%   Y : DIC计算结果的X坐标
%   U : DIC计算的水平位移
%   V : DIC计算的竖直位移
%   ZNCC : DIC计算结果的相关系数
%   S : DIC计算结果的有效位

prefix = ['Polynomial/L', num2str(L), 'N', num2str(N), 'M', num2str(M)];

x = readmatrix([prefix, '_x.csv']);
y = readmatrix([prefix, '_y.csv']);
u = readmatrix([prefix, '_u.csv']);
v = readmatrix([prefix, '_v.csv']);
zncc = readmatrix([prefix, '_zncc.csv']);
s = readmatrix([prefix, '_valid_sign.csv']);

end