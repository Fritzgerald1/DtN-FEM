[Ko,Kl,Kr,l,r,~] = Stiffness_Changed_Position(Kg,bnd_L,bnd_R,Nnode); % 交换刚度矩阵位置
%% 左边界入射模态的节点位移/力
bnd = bnd_L;
n = mode_in;
[P_inc,~] = nodal_displacement_and_force( wd,kd,n,Amp,bnd,Coordinate,lambda,mu,density,h,0 );
% P_inc = -P_inc;
% P_inc(1:2:end) = -P_inc(1:2:end);
% P_inc(2:2:end) = -P_inc(2:2:end);
%% 右边界透射模态的节点位移/力
bnd = bnd_R;
n = mode_in;
[P_R,~] = nodal_displacement_and_force( wd,kd,n,Amp,bnd,Coordinate,lambda,mu,density,h,0 );
% P_R = -P_R;
% P_R(1:2:end) = -P_R(1:2:end);
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
u_Tx = u_T(1:2:end);
u_Ty = u_T(2:2:end);

figure(Position=[263,382,1106,420])
subplot(2,1,1)
px = Coordinate(bnd_T,1);
[px1,I] = sort(px);
u_Tx = u_Tx(I);
u_Ty = u_Ty(I);
plot(px1,real(u_Tx),'b-',LineWidth=1.5)
hold on
plot(px1,imag(u_Tx),'b-.')
hold off
title('u')
legend('real','image')
subplot(2,1,2)
plot(px1,real(u_Ty),'r-',LineWidth=1.5)
hold on
plot(px1,imag(u_Ty),'r-.')
hold off
title('v')
legend('real','image')


% plotWave(1.5e-3,Coordinate,u)


function plotWave(x,Coordinate,u)
X = Coordinate(:,1);
bnd = find(X==x);
uB = u(bnd);
U = uB(1:2:end);
V = uB(2:2:end);
py = Coordinate(bnd,2);

plot(real(U),py,'b-',LineWidth=1.5);
hold on
plot(imag(U),py,'b-.')
hold off
end
