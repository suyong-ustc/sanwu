function [x, y, u, v, zncc, s] = ReadSinDICData(T, M, U, N)
%ReadSinDICData 读取正弦变形的DIC分析结果
%   T : 正弦变形周期
%   M : 子区尺寸
%   U : 空缺
%   N : 形函数阶数
%   X : DIC计算结果的X坐标
%   Y : DIC计算结果的X坐标
%   U : DIC计算的水平位移
%   V : DIC计算的竖直位移
%   ZNCC : DIC计算结果的相关系数
%   S : DIC计算结果的有效位

prefix = ['Sinusoidal/T', num2str(T), 'M', num2str(M), 'U', num2str(U), 'N', num2str(N)];

x = readmatrix([prefix, '_x.csv']);
y = readmatrix([prefix, '_y.csv']);
u = readmatrix([prefix, '_u.csv']);
v = readmatrix([prefix, '_v.csv']);
zncc = readmatrix([prefix, '_zncc.csv']);
s = readmatrix([prefix, '_valid_sign.csv']);

end