[Ko,Kl,Kr,l,r,~] = Stiffness_Changed_Position(Kg,bnd_L,bnd_R,Nnode); % 交换刚度矩阵位置
%% 左边界入射模态的节点位移/力
bnd = bnd_L;
n = mode_in;
[P_inc,~] = nodal_displacement_and_force( wd,kd,n,Amp,bnd,Coordinate,lambda,mu,density,h,0 );
% P_inc = conj(P_inc);
% P_inc = -P_inc;
% P_inc(1:2:end) = -P_inc(1:2:end);
%% 右边界透射模态的节点位移/力
bnd = bnd_R;
n = mode_in;
[P_R,~] = nodal_displacement_and_force( wd,kd,n,Amp,bnd,Coordinate,lambda,mu,density,h,180 );
% P_R = conj(P_R);
% P_R = -P_R;

% P_R(2:2:end) = -P_R(2:2:end);
%% 等效动力刚度矩阵
KPT = sparse(Ko);
%% 等效节点力
c_inc = 1; % 令入射系数c_inc为1
c_TR = 1; % 令右端透射系数为1
F = -( Kl*P_inc*c_inc + Kr*P_R*c_TR);
% F = -(Kr*P_R*c_TR);
%% 求解 KPT*[u;coeff]=F
sol = KPT\sparse(F);
%% 判断矩阵KPT是否病态
% cond(KPT) % KPT的条件数
% dKPT = decomposition(KPT);
% isIllConditioned(dKPT)
%% 左右端位移
u = nan(2*Nnode,1);
u(l) = P_inc*c_inc;
u(r) = P_R*c_TR;
u(isnan(u)) = sol;
%% 绘画出上表面的位移
bottom = sort([2*bnd_B-1;2*bnd_B]);
u_T = u(bottom);
uT = u_T(1:2:end);
vT = u_T(2:2:end);
figure(Position=[263,382,1106,420])
subplot(2,1,1)
px = Coordinate(bnd_T,1);
[px1,I] = sort(px);
uT = uT(I);
vT = vT(I);
plot(px1,real(uT),'b-',LineWidth=1.5)
hold on
plot(px1,imag(uT),'b-.')
hold off
title('u')
legend('real','image')
subplot(2,1,2)
plot(px1,real(vT),'r-',LineWidth=1.5)
hold on
plot(px1,imag(vT),'r-.')
hold off
title('v')
legend('real','image')


plotWave(0,Coordinate,u)
plotWave(15e-3,Coordinate,u)

% function plotWave(x,Coordinate,u)
% X = Coordinate(:,1);
% bnd_node = find(X==x);
% bnd = sort([2*bnd_node-1;2*bnd_node]);
% uB = u(bnd);
% U = uB(1:2:end);
% V = uB(2:2:end);
% py = Coordinate(bnd_node,2);
% [py1,I] = sort(py);
% U = U(I);
% V = V(I);
% 
% figure
% subplot(1,2,1)
% plot(real(U),py1,'b-',LineWidth=1.5);
% hold on
% plot(imag(U),py1,'b-.')
% hold off
% subtitle('u1')
% 
% subplot(1,2,2)
% plot(real(V),py1,'b-',LineWidth=1.5);
% hold on
% plot(imag(V),py1,'b-.')
% hold off
% subtitle('u2')
% end
