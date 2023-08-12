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
	[P_Rn,T_Rn] = nodal_displacement_and_force( wd,kd,n,Amp,bnd,Coordinate,lambda,mu,density,h,pi/2 );
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
[r,c] = size(T1);
T(end-r+1:end,end-c+1:end) = T1;

KPT = sparse(KP+T);
% KPi = Kl*P_inc;
% Ti = zeros(size(KPi));
% Ti = zeros(2*Nnode,1);
Ti = sparse(2*Nnode,1);
[r,c] = size(T_inc);
Ti(end-r+1:end,end-c+1:end) = T_inc;
F = sparse(-(Kl*P_inc+Ti)); % 令入射cL为1e-6

%% 求解 KPT*[u;coeff]=F
sol = KPT\sparse(F);

%% 散射系数
coeff = sol(end-2*modes+1:end);
c_RE = coeff(1:modes);
c_TR = coeff(modes+1:end);

cond(KPT) % KPT的条件数
dKPT = decomposition(KPT);
isIllConditioned(dKPT) % 是否是病态矩阵


% coeff = Us(end-c+1:end,:);