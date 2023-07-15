function [V] = get_amp(kd,wd)
%% 材料属性
CL = 6.35;
CT = 3.13;
h=0.5;
mu = 26;

%% 量纲还原
k = kd/h;
w = wd*CT/h;
p = sqrt(w.^2/CL^2 - k.^2);
q = sqrt(w.^2/CT^2 - k.^2);

%% 频散方程
M(1,1) = mu*(-2i*k*p*sin(p*h));
M(1,2) = mu*((k^2-q^2)*sin(q*h));
M(2,1) = -mu*((q^2-k^2)*cos(p*h));
M(2,2) = -2*mu*1i*k*q*cos(q*h);
[Val,~] = eig(M);
V = Val(:,1);
% 
% [L,U] = lu(M);
% z = [1;0];
% V = U\(L\z);

end