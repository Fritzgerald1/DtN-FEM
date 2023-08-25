clear
close all
%% 子文件夹路径
addpath('.\import'); % 读入网格信息
addpath('.\K'); % 形成刚度矩阵
addpath('.\lamb') % Lamb波波数、频率及幅值
%% Aluminum材料参数
density = 2.7*10^3;
lambda = 51e9;
mu = 26e9;
% CL = 6.35e3; % 纵波波速
% CT = 3.13e3; % 横波波速
CL = sqrt((lambda+2*mu)/density);
CT = sqrt(mu/density);
d = 1e-3; % 板厚
h = d/2;
ff = 1e6:.1e6:2e6;
nff = length(ff);
cR = zeros(2,nff);
cT = zeros(2,nff);
%% 几何网格
% Import Mesh Information %
% bnd_401是左边界节点编号，bnd_501是右边界，bnd_101是自由应力边界（上下边界及缺陷边界） %
[Nnode, Nelement, Coordinate, Ielement] = read_mesh_info("mesh.dat");
[bnd_L, bnd_R, bnd_T,bnd_B,bnd_free, bnd_L_e, bnd_R_e, bnd_free_e] = bound_information(Coordinate, Ielement);
%% 刚度矩阵
[K,M] = Stiffness_Mass_matrix( Nnode,Nelement,Ielement,Coordinate,lambda,mu,density);
%% DtN-FEM
for tt = 1:nff
	%% 入射波
	f = ff(tt); % 入射频率
	w = 2*pi*f;
	mode_in = 1; % A0模态入射
	%% 动力刚度矩阵
	Kg = sparse((K-w^2*M));
	%% Lamb波频散特性
	wd = w/CT*h; % 无量纲频率
	Fun = @lamb;
	[kd,~,modes] = get_wavenumber(wd,lambda,mu,density,h,Fun); % 得到无量纲波数和模态个数
	% kd = kd(end:-1:1); % 波数按低阶到高阶排列
	k = kd/h; % 去归一化
	[Amp] = get_amplitude( kd,wd,lambda,mu,density,h,Fun );
	%% 得到散射系数和位移
	% c_RE:	反射系数；
	% c_TR:	透射系数；
	% u:	位移;
	[c_RE,c_TR,u] = PT(bnd_L,bnd_R,mode_in,modes,wd,kd,Amp,Coordinate,Nnode,lambda,mu,density,h, Kg,bnd_T);
	cR(1:2,tt) = full(c_RE(1:2));
	cT(1:2,tt) = full(c_TR(1:2));
end

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