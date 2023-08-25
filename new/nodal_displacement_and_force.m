function [U,T] = nodal_displacement_and_force(wd,Kd,n,Amp,bnd,Coordinate,lambda,mu,density,h, alpha)
%% 左右边界处的节点模态节点位移Φ和力T (第n阶模态)
CL = sqrt((lambda+2*mu)/density);
CT = sqrt(mu/density); 

kd = Kd(n); % 第k阶波数
k = kd/h;
w = wd*CT/h;
f = w/2/pi;

p = sqrt((w/CL)^2-k.^2);
q = sqrt((w/CT)^2-k.^2);

x = Coordinate(bnd,1);
y = Coordinate(bnd,2);
dy = y(2)-y(1);

E = [exp(1i*p*y) exp(-1i*p*y) exp(1i*q*y) exp(-1i*q*y)];
A = Amp(:,n); % 第n阶模态的幅值系数

u =  E * diag([k, k, -q, q].*1i) * A .* exp(1i*k*x);
v =  E * diag([p, -p k k].*1i) * A .* exp(1i*k*x);
% u = imag(u); v = real(v);

Ux = cos(alpha)*u + sin(alpha)*v;
Uy = sin(alpha)*u + cos(alpha)*v;

t11 = E * diag([mu*(2*p^2-q^2-k^2)*[1 1], 2*mu*q*k*[1 -1]]) * A .* exp(1i*k*x) .*dy;
t12 = E * diag([2*mu*k*p*[-1 1], mu*(q^2-k^2)*[1 1]]) * A .* exp(1i*k*x) .*dy;
t22 = E * diag([-mu*(q^2-k^2)*[1 1], 2*mu*k*q*[-1 1]]) * A .* exp(1i*k*x) .*dy;
% t11 = real(t11); t12 = imag(t12); t22 = real(t22); 
tx = t11+t12;
ty = t22+t12;

% Mx = diag([1i*mu*k*(2*p^2-q^2-k^2-2*k*p), 1i*mu*k*(2*p^2-q^2-k^2+2*k*p), 1i*mu*k*(q^2-k^2+2*k*q), 1i*mu*k*(q^2-k^2-2*k*q)]);
% My = diag([-1i*mu*p*(2*k*p+(q^2-k^2)), -1i*mu*p*(2*k*p-(q^2-k^2)), -1i*mu*q*(2*q*k-(q^2-k^2)),  1i*mu*q*(2*q*k+(q^2-k^2))]);
% tx = E * Mx * A .* exp(1i*k*x);
% ty = E * My * A .* exp(1i*k*x);
% tx = real(tx); ty = real(ty);

Tx = cos(alpha)*tx + sin(alpha)*ty;
Ty = -sin(alpha)*tx + cos(alpha)*ty;

% 合并x和y的位移/力
tmpU = [Ux,Uy]';
U = tmpU(:);

tmpT = [Tx,Ty]';
T = tmpT(:);

%% 绘图检验
% 位移
% figure
% subplot(1,2,1)
% [y1,I] = sort(y);
% plot(u(I),y1)
% grid on
% hold on
% plot(v(I),y1,'--')
% legend('u','v')
% hold off
% 
% % 外力
% subplot(1,2,2)
% [y1,I] = sort(y);
% plot(tx(I),y1)
% grid on
% hold on
% plot(ty(I),y1,'--')
% legend('tx','ty')
% sgtitle(['第' num2str(n) '阶模态'])
% hold off
% 
% % 应力
% figure
% [y1,I] = sort(y);
% plot(t11(I),y1)
% grid on
% hold on
% plot(t12(I),y1,'--')
% plot(t22(I),y1,'.-')
% legend('t11','t12','t22')
% sgtitle(['第' num2str(n) '阶模态'])
% hold off
end