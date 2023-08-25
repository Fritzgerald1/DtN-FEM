function [c_RE,c_TR,u] = PT(bnd_L,bnd_R,mode_in,modes,wd,kd,Amp,Coordinate,Nnode,lambda,mu,density,h, Kg,bnd_T)

[Ko,Kl,Kr,l,r,~] = Stiffness_Changed_Position(Kg,bnd_L,bnd_R,Nnode); % 交换刚度矩阵位置
%% 左边界入射模态的节点位移/力
bnd = bnd_L;
n = mode_in;
[P_inc,T_inc] = nodal_displacement_and_force( wd,kd,n,Amp,bnd,Coordinate,lambda,mu,density,h,0 );


%% 左边界反射模态的节点位移/力
bnd = bnd_L;
P_L = zeros(2*length(bnd_L),modes);
T_L = zeros(2*length(bnd_L),modes);
for n = 1:modes	
	[P_Ln,T_Ln] = nodal_displacement_and_force( wd,kd,n,Amp,bnd,Coordinate,lambda,mu,density,h,0 );
	P_L(:,n) = P_Ln;
	T_L(:,n) = T_Ln;
end

%% 右边界透射模态的节点位移/力
bnd = bnd_R;
P_R = zeros(2*length(bnd_R),modes);
T_R = zeros(2*length(bnd_R),modes);

for n = 1:modes
	[P_Rn,T_Rn] = nodal_displacement_and_force( wd,kd,n,Amp,bnd,Coordinate,lambda,mu,density,h,0 );
	P_R(:,n) = P_Rn;
	T_R(:,n) = T_Rn;
end

%% 用系数代替边界节点位移
KP = [Ko, Kl*P_L, Kr*P_R];

T = zeros(size(KP));
% T1 = [T_L, zeros(size(T_R));
% 	zeros(size(T_L)), T_R] ;
T1 = blkdiag(T_L,T_R);
% r = 2*(bnd_L+bnd_R); c = 1+nmode*2;
[R,C] = size(T1);
T(end-R+1:end,end-C+1:end) = T1;
T = sparse(T);

KPT = KP+T;
Ti = sparse(2*Nnode,1);
[R,C] = size(T_inc);
Ti(end-R+1:end,end-C+1:end) = T_inc;

c_inc = 1e-6; % 令入射系数c_inc为1e-6
F = sparse(-(Kl*P_inc+Ti)*c_inc);

%% 判断矩阵KPT是否病态
% cond(KPT) % KPT的条件数
% dKPT = decomposition(KPT);
% isIllConditioned(dKPT)

%% 求解 KPT*[u;coeff]=F
sol = KPT\sparse(F);
coeff = sol(end-2*modes+1:end);
c_RE = coeff(1:modes)./c_inc;
c_TR = coeff(modes+1:end)./c_inc;

u = zeros(2*Nnode,1);
u(l)=P_inc*c_inc+P_L*c_RE;
u(r)=P_R*c_TR;

u(u==0) = sol(1:end-2*modes);

%% 能量检验
uI = P_inc*c_inc;
ul = P_L*c_RE;
ur = P_R*c_TR;
tI = T_inc*c_inc;
tl = T_L*c_RE;
tr = T_R*c_TR;

eIN = sum(abs(uI.*tI));
eOUT = sum(ur.*tr)+sum(ur.*tr);

% energyIn = uI
% energyOut = 
%% 绘画出上表面的位移
% top = sort([2*bnd_T-1;2*bnd_T]);
% u_T = u(top);
% u_Tx = u_T(1:2:end);
% u_Ty = u_T(2:2:end);
% figure
% subplot(2,1,1)
% plot(linspace(0,15e-3,length(bnd_T)),u_Tx,'b')
% title('u')
% subplot(2,1,2)
% plot(linspace(0,15e-3,length(bnd_T)),u_Ty,'r')
% title('v')

end
