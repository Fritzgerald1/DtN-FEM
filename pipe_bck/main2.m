clear
%% 子文件夹路径
addpath(genpath('.\import')); % 读入网格信息
addpath(genpath('.\K')); % 形成刚度矩阵
addpath(genpath('.\lamb')) % Lamb波波数、频率及幅值

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

%% 入射波
f = 1e6; % 入射频率
w = 2*pi*f;
mode_in = 1; % A0模态入射

%% 几何网格
% Import Mesh Information %
% bnd_401是左边界节点编号，bnd_501是右边界，bnd_101是自由应力边界（上下边界及缺陷边界） %
[Nnode, Nelement, Coordinate, Ielement] = read_mesh_info("mesh_pipe.dat");
[bnd_L, bnd_R,bnd_T,bnd_B,bnd_free, bnd_L_e, bnd_R_e, bnd_free_e] = bound_information(Coordinate, Ielement);

%% 动力刚度矩阵
% K: global stiffness matrix %
% M: global mass matrix %
% Kg: Dynamic Stiffness Matrix
[Kg,K] = Stiffness_Mass_matrix( Nnode,Nelement,Ielement,Coordinate,lambda,mu,density,w);
[Ko,Kl,Kr] = Stiffness_Changed_Position(Kg,bnd_L,bnd_R,Nnode); % 交换刚度矩阵位置

%% Lamb波频散特性
wd = w/CT*h; % 无量纲频率
Fun = @lamb;
[kd,~,modes] = get_wavenumber(wd,lambda,mu,density,h,Fun); % 得到无量纲波数和模态个数
kd = kd(end:-1:1); % 波数按低阶到高阶排列
k = kd/h; % 去归一化
[Amp] = get_amplitude( kd,wd,lambda,mu,density,h,Fun );
