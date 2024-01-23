function [xx, theory_u] = PolyDisplacementByTheory(x, L, N, M)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
    xx = min(x) : 0.1 : max(x)+0.5;
    
    M = (M - 1)/2;
    zeta = xx / (2*M+1);
    
    u0 = 2^-(L+1) .* (1 - (2*zeta-1).^(L+1)) ./ ((L+1)*(1-zeta));
    
    u1 = u0 + 3*2^-(L+1) .* zeta./(1-zeta) .* ...
        ( (1+(2*zeta-1).^(L+1)) ./ ((L+1)*(1-zeta)) - (1-(2*zeta-1).^(L+2)) ./ ((L+1)*(L+2)*(1-zeta).^2) );
    
    u2 = u1 + 5*2^-(L+1) * (1.5*(zeta./(1-zeta)).^2 - 0.5) .* ...
        ( (1-(2*zeta-1).^(L+1)) ./ ((L+1)*(1-zeta)) - ...
          3 * (1+(2*zeta-1).^(L+2)) ./ ((L+1)*(L+2)*(1-zeta).^2) + ...
          3 * (1-(2*zeta-1).^(L+3)) ./ ((L+1)*(L+2)*(L+3)*(1-zeta).^3));
    
    ratio = 0;
    if L == 2
        a = 1e-3;
        ratio = a * (2*M+1)^2;
    elseif L == 3
        a = 1e-4;
        ratio = a * (2*M+1)^3;
    elseif L == 4
        a = 1e-6;
        ratio = a * (2*M+1)^4;
    elseif L == 5
        a = 1e-7;
        ratio = a * (2*M+1)^5;
    end
    
    if N == 0
        theory_u = u0 * ratio;
    elseif N == 1
        theory_u = u1 * ratio;
    else
        theory_u = u2 * ratio;
    end

end