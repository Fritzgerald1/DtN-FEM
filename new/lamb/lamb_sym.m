function [error,varargout] = lamb_sym(kd,wd,lambda,mu,density,h)
%% 材料属性
CL = sqrt((lambda+2*mu)/density);
CT = sqrt(mu/density); 

%% 量纲还原
k = kd/h;
w = wd*CT/h;
p = sqrt(w.^2/CL^2 - k.^2);
q = sqrt(w.^2/CT^2 - k.^2);

%% 频散方程
M(1,1) = 2i*k*p*sin(p*h);
M(1,2) = (k^2-q^2)*sin(q*h);
M(2,1) = (q^2-k^2)*cos(p*h);
M(2,2) = -2i*k*q*cos(q*h);
error = det(M);

if nargout == 2
	flag=1;
	if abs(p)<1e-1,	flag=0;	end
	if abs(q)<1e-1,	flag=0;	end
	varargout{1} = flag;
end

if nargout == 3
	flag=1;
	varargout{1}=flag;
	
	V = null(M);
	V = V(:,1)/V(1,1);
	varargout{2} = V;
end
