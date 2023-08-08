function [error,varargout] = lamb(kd,wd)
%% 材料属性
CL = 6.35;
CT = 3.13;
h = 0.5;
mu = 26;
lambda = 51;

%% 量纲还原
k = kd/h;
w = wd*CT/h;
KL = sqrt(w.^2/CL^2 - k.^2);
KT = sqrt(w.^2/CT^2 - k.^2);

%% 频散方程
M(1,1) = -2*mu*k*KL*exp(1i*KL*h);
M(1,2) = 2*mu*k*KL*exp(-1i*KL*h);
M(1,3) = mu*(KT^2-k^2)*exp(1i*KT*h);
M(1,4) = mu*(KT^2-k^2)*exp(-1i*KT*h);
M(2,1) = -mu*(KT^2-k^2)*exp(1i*KL*h);
M(2,2) = -mu*(KT^2-k^2)*exp(-1i*KL*h);
% M(2,1) = (-lambda*(w/CL)^2-2*mu*KL^2)*exp(1i*KL*h);
% M(2,2) = (-lambda*(w/CL)^2-2*mu*KL^2)*exp(-1i*KL*h);
M(2,3) = -2*mu*k*KT*exp(1i*KT*h);
M(2,4) = 2*mu*k*KT*exp(-1i*KT*h);

M(3,1) = -2*mu*k*KL*exp(-1i*KL*h);
M(3,2) = 2*mu*k*KL*exp(1i*KL*h);
M(3,3) = mu*(KT^2-k^2)*exp(-1i*KT*h);
M(3,4) = mu*(KT^2-k^2)*exp(1i*KT*h);
M(4,1) = -mu*(KT^2-k^2)*exp(-1i*KL*h);
M(4,2) = -mu*(KT^2-k^2)*exp(1i*KL*h);
% M(4,1) = (-lambda*(w/CL)^2-2*mu*p^2)*exp(-1i*p*h);
% M(4,2) = (-lambda*(w/CL)^2-2*mu*p^2)*exp(1i*p*h);
M(4,3) = -2*mu*k*KT*exp(-1i*KT*h);
M(4,4) = 2*mu*k*KT*exp(1i*KT*h);

error = det(M);

if nargout == 2
	flag=1;
	if abs(KL)<1e-1,	flag=0;	end
	if abs(KT)<1e-1,	flag=0;	end
	varargout{1} = flag;
end

if nargout == 3
	flag=1;
	varargout{1}=flag;
	
	V = null(M,1);
	V = V(:,1)/V(1,1);
	varargout{2} = V;
end
