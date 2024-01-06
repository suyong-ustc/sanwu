function ShowSinDICData(T, M, U, N)
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

[x, y, u, v, zncc, s] = ReadSinDICData(T, M, U, N);

xx = x(1,:);
yy = y(:,1)';

figure;
subplot(2,2,1);
imagesc(xx,yy,u);
axis equal;
colorbar;

subplot(2,2,2);
imagesc(xx,yy,v);
axis equal;
colorbar;

subplot(2,2,3);
imagesc(xx,yy,zncc);
axis equal;
colorbar;

subplot(2,2,4);
imagesc(xx,yy,s);
axis equal;
colorbar;

end