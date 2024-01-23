clc;
clear;
close all;


%% 分析不同阶的变形

L = 3;
M = 49;

[~, u0_mean, u0_std] = PolyDisplacementByDIC(L, 0, M);
[~, u1_mean, u1_std] = PolyDisplacementByDIC(L, 1, M);
[x, u2_mean, u2_std] = PolyDisplacementByDIC(L, 2, M);

[~, theory_u0] = PolyDisplacementByTheory(x, L, 0, M);
[~, theory_u1] = PolyDisplacementByTheory(x, L, 1, M);
[xx, theory_u2] = PolyDisplacementByTheory(x, L, 2, M);

if L == 2
    real_u = 1e-3 * xx.^2;
elseif L == 3
    real_u = 1e-4 * xx.^3;
elseif L == 4
    real_u = 1e-6 * xx.^4;
else
    real_u = 1e-7 * xx.^5;
end

figure;
hold on;
errorbar(x, u0_mean, u0_std,'r');
errorbar(x, u1_mean, u1_std,'b');
errorbar(x, u2_mean, u2_std,'g');
plot(xx, theory_u0,'r');
plot(xx, theory_u1,'b');
plot(xx, theory_u2,'g');
plot(xx, real_u,'k');
xlim([0, (M+1)/2]);

data1 = [x', u0_mean', u0_std', u1_mean', u1_std', u2_mean', u2_std'];
data2 = [xx', theory_u0', theory_u1', theory_u2', real_u'];