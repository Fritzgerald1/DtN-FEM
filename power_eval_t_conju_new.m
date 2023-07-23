function [ Pmn ] = power_eval_t_conju_new( imode,n_layer,mu,wd,kd,zamp_mode,x1,ipm )
% single layer
% 输入变量：
	% imode - 第i阶
	% wd - 无量纲频率
	% kih - 第i阶无量纲波数
y = [-1 0 1];

%ccd
n = 8;
[x,w]=gauss_integration(n);
for i =1:n_layer
    rmui=mu(1,i)/mu(1,1); 
    y2=y(i+1);
    y1=y(i);
    he = 0.5*(y2-y1);
    p = sqrt(wd^2-kd^2);
	%% 高斯积分
    for iw=1:n
        ys1=x1;
        ys2=((y2-y1)*x(iw)+y1+y2)*0.5;
        z_um1=(zamp_mode(2*i-1,imode)*exp(1i*p*ys2)+zamp_mode(2*i,imode)*exp(-1i*p*ys2))*exp(ipm*1i*kd*ys1);
        z_um2=(zamp_mode(2*i-1,imode)*exp(1i*p*ys2)+zamp_mode(2*i,imode)*exp(-1i*p*ys2))*exp(-ipm*1i*kd*ys1);
        z_um2=conj(z_um2);
        Pmn= Pmn-real(0.5*kd*rmui*z_um1*z_um2*he*w(iw));
    end
end
%Pmn
end