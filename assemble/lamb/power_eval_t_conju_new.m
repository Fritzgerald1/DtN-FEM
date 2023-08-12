function [ Pm ] = power_eval_t_conju_new( imode,wd,kd,Amp)
% imode - 第i阶
% wd - 无量纲频率
% kih - 第i阶无量纲波数
n = 8;
[x,w]=gauss_integration(n);
he = 1;
p = sqrt((wd*CT/CL)^2-kd^2);
q = sqrt(wd^2-kd^2);
A = Amp(:,imode);
%% 高斯积分
for iw=1:n
	x2 = x(iw);
	[u1,u2,t11,t12,t22] = ut_lamb(x2);
	% 	z_um1=(A*exp(1i*p*x2)+Amp(2,imode)*exp(-1i*p*x2))*exp(ipm*1i*kd*x1); % 一层的z位移
	% 	z_um2=(Amp(1,imode)*exp(1i*p*x2)+Amp(2,imode)*exp(-1i*p*x2))*exp(-ipm*1i*kd*x1); % 二层的z位移
	% 	z_um2=conj(z_um2);

	Pm = Pm + real(); % 单位能量
	Pm= Pm-real(0.5*kd*z_um1*z_um2*he*w(iw)); % P = int{ σ13 * (u3)', dx2 } = int{ 1/2 * kh * (uT*u) , dx2}
end
%Pmn
end

function [u1,u2,t11,t12,t22] = ut_lamb(x2)
E = [exp(1i*p*x2) exp(-1i*p*x2) exp(1i*q*x2) exp(-1i*q*x2)];
u1 =  E * diag([k, k, -q, q].*1i) * A ;
u2 =  E * diag([p, -p k k].*1i) * A ;
t11 = E * diag([mu*(2*p^2-q^2-k^2)*[1 1], 2*mu*q*k*[1 -1]]) * A;
t12 = E * diag([2*mu*k*p*[-1 1], mu*(q^2-k^2)*[1 1]]) * A;
t22 = E * diag([-mu*(q^2-k^2)*[1 1], -2*mu*k*q*[-1 1]]) * A;
end