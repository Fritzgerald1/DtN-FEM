function [U,T] = nodal_displacement_and_force(wd,Kd,n,Amp,bnd,Coordinate,lambda,mu,density,h, alpha)
%% ���ұ߽紦�Ľڵ�ģ̬�ڵ�λ�Ʀ�����T (��n��ģ̬)
CL = sqrt((lambda+2*mu)/density);
CT = sqrt(mu/density); 

kd = Kd(n); % ��k�ײ���
k = kd/h;
w = wd*CT/h;

p = sqrt((w/CL)^2-k.^2);
q = sqrt((w/CT)^2-k.^2);

x = Coordinate(bnd(1),1);
% x = linspace(0,15e-3,1000);
y = Coordinate(bnd,2);
% y = y(1);


E = [exp(1i*p*y) exp(-1i*p*y) exp(1i*q*y) exp(-1i*q*y)];
A = Amp(:,n); % ��n��ģ̬�ķ�ֵϵ��

u =  E * diag([k, k, -q, q].*1i) * A .* exp(1i*k*x);
v =  E * diag([p, -p k k].*1i) * A .* exp(1i*k*x);
u = imag(u); v = real(v);

Ux = cos(alpha)*u + sin(alpha)*v;
Uy = sin(alpha)*u + cos(alpha)*v;
% u = u./max(abs(u)); v = v./max(abs(v));

t11 = E * diag([mu*(2*p^2-q^2-k^2)*[1 1], 2*mu*q*k*[1 -1]]) * A .* exp(1i*k*x);
t12 = E * diag([2*mu*k*p*[-1 1], mu*(q^2-k^2)*[1 1]]) * A .* exp(1i*k*x);
t22 = E * diag([-mu*(q^2-k^2)*[1 1], 2*mu*k*q*[-1 1]]) * A .* exp(1i*k*x);
t11 = real(t11); t12 = imag(t12); t22 = real(t22); 
tx = t11+t12;
ty = t12+t22;
Tx = cos(alpha)*tx + sin(alpha)*ty;
Ty = sin(alpha)*tx + cos(alpha)*ty;

%% �ϲ�x��y��λ��/��
tmpU = [Ux,Uy]';
U = tmpU(:);

tmpT = [Tx,Ty]';
T = tmpT(:);



% %% ��ͼ����
% %λ��
% figure;
% [y1,I] = sort(y);
% plot(u(I),y1)
% grid on
% hold on
% plot(v(I),y1,'--')
% legend('u','v')
% %��
% figure;
% [y1,I] = sort(y);
% plot(tx(I),y1)
% grid on
% hold on
% plot(ty(I),y1,'--')
% legend('tx','ty')

end