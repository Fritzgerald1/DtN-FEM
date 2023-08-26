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
% u = -imag(u); v = real(v);

u = cosd(alpha)*u + sind(alpha)*v;
v = -sind(alpha)*u + cosd(alpha)*v;

% t11 = E * diag([(lambda*(-p^2-k^2)-2*mu*k^2)*[1 1], 2*mu*q*k*[1 -1]]) * A .* exp(1i*k*x) .*dy;
% t12 = E * diag([2*mu*k*p*[-1 1], mu*(q^2-k^2)*[1 1]]) * A .* exp(1i*k*x) .*dy;
% t22 = E * diag([(lambda*(-p^2-k^2)-2*mu*p^2)*[1 1], 2*mu*k*q*[-1 1]]) * A .* exp(1i*k*x) .*dy;

t11 = E * diag([mu*(2*p^2-q^2-k^2)*[1 1], 2*mu*q*k*[1 -1]]) * A .* exp(1i*k*x) .*dy;
t12 = E * diag([2*mu*k*p*[-1 1], mu*(q^2-k^2)*[1 1]]) * A .* exp(1i*k*x) .*dy;
t22 = E * diag([-mu*(q^2-k^2)*[1 1], 2*mu*k*q*[-1 1]]) * A .* exp(1i*k*x) .*dy;
% t11 = real(t11); t12 = -imag(t12); t22 = real(t22); 
tx = t11;
ty = t12;

tx = cosd(alpha)*tx + sind(alpha)*ty;
ty = -sind(alpha)*tx + cosd(alpha)*ty;

tmpU = [u,v]'; U = tmpU(:); % 合并x和y的位移
tmpT = [tx,ty]'; T = tmpT(:); % 合并x和y的力
%% 绘图检验
%{

figure(Position=[263,382,1106,420])
[yi,I] = sort(y);
y1 = yi./max(yi);
% 位移

ax1 = subplot(1,2,1);
plot(u(I)./max(abs([u;v])),y1,'LineWidth',1)
axis([-1.1 1.1 -1.1 1.1])
grid on
hold on
plot(v(I)./max(abs([u;v])),y1,'--','LineWidth',1)
legend('u','v')
ax1.XAxisLocation='origin';ax1.YAxisLocation='origin';
hold off

% 外力
ax2 = subplot(1,2,2);
plot(tx(I)./max(abs([tx;ty])),y1,'LineWidth',1)
axis([-1.1 1.1 -1.1 1.1])
grid on
hold on
plot(ty(I)./max(abs([tx;ty])),y1,'--','LineWidth',1)
legend('tx','ty')
ax2.XAxisLocation='origin';ax2.YAxisLocation='origin';

str_title = sprintf('第%d阶模态\n k = %.0fm^{-1},  f = %.2fMHz',n,k,f/1e6);
sgtitle(str_title)
hold off

% 应力
figure(Position=[263,382,1106,420])
[y1,I] = sort(y);
plot(t11(I),y1,'LineWidth',1)
grid on
hold on
plot(t12(I),y1,'--','LineWidth',1)
plot(t22(I),y1,'.-','LineWidth',1)
legend('t11','t12','t22')
ax = gca;
ax.XAxisLocation='origin';ax.YAxisLocation='origin';
str_title = sprintf('第%d阶模态\n k = %.0fm^{-1},  f = %.2fMHz',n,k,f/1e6);
title(str_title)
hold off

%}
end