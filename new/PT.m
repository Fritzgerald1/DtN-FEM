[Ko,Kl,Kr,l,r,~] = Stiffness_Changed_Position(Kg,bnd_L,bnd_R,Nnode); % 交换刚度矩阵位置
%% 左边界入射模态的节点位移/力
bnd = bnd_L;
n = mode_in;
[P_inc,~] = nodal_displacement_and_force( wd,kd,n,Amp,bnd,Coordinate,lambda,mu,density,h,0 );
%% 右边界透射模态的节点位移/力
bnd = bnd_R;
n = mode_in;
[P_R,~] = nodal_displacement_and_force( wd,kd,n,Amp,bnd,Coordinate,lambda,mu,density,h,0 );
% P_R = -P_R;
%% 等效动力刚度矩阵
KPT = sparse(Ko);
%% 等效节点力
F = Kl*P_inc + Kr*P_R;
% c_inc = 1; % 令入射系数c_inc为1
% c_TR = 1; % 令右端透射系数为1

%% 判断矩阵KPT是否病态
% cond(KPT) % KPT的条件数
% dKPT = decomposition(KPT);
% isIllConditioned(dKPT)

%% 求解 KPT*[u;coeff]=F
sol = KPT\sparse(F);
u = nan(2*Nnode,1);

u(l) = P_inc*c_inc;
u(r) = P_R*c_TR;
u(isnan(u)) = sol;

%% 绘画出上表面的位移
top = sort([2*bnd_T-1;2*bnd_T]);
u_T = u(top);
u_Tx = u_T(1:2:end);
u_Ty = u_T(2:2:end);

figure(Position=[263,382,1106,420])
subplot(2,1,1)
plot(linspace(0,max(Coordinate(:,1)),length(bnd_T)),u_Tx,'b')
title('u')
subplot(2,1,2)
plot(linspace(0,max(Coordinate(:,1)),length(bnd_T)),u_Ty,'r')
title('v')

% end
