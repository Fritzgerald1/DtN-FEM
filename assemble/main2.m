% FEM SH wave %
% Input %
% density:matterial parameter %
% h:unit length of x3 direction %
% mu:shear modulus %
% Cf:cicular frequency %

%% 子文件夹路径
addpath(genpath('.\import')); % 读入网格信息
addpath(genpath('.\K')); % 形成刚度矩阵
addpath(genpath('.\lamb')) % Lamb波波数、频率及幅值

%% Aluminum材料参数
density = 2.7*10^3; 
mu = 26.1*10^9; 
CT = 3.13e3; % 剪切波速
% CT = sqrt(mu1/density1); 
d = 1e-3; % 板厚
h = d/2;

%% A0模态入射
f = 1.5*10^6; % 入射频率
cf = 2*pi*f;
mode_in = 2; 

%% Geometry Setting
L=1.5*10^(-2);
W=1*10^(-3);
L_crack = W;

% Import Mesh Information %
% bnd_401是左边界节点编号，bnd_501是右边界，bnd_101是自由应力边界（上下边界及缺陷边界） %
[Nnode, Nelement, Coordinate, Ielement, dL, dW] = read_mesh_info();
[bnd_401, bnd_501, bnd_101, bnd_401_e, bnd_501_e, bnd_101_e] = bound_information(Coordinate, Ielement,L,W,L_crack);

Nelement_crack = round(L_crack/dL);
Nnode_crack = 2*Nelement_crack-1;
L_ele = round(L/dL);
L_nod_s = 2*L_ele+1;
L_nod_d = L_ele+1;

%% Stiffness matrix
% K: global stiffness matrix %
% M: global mass matrix %
% Kg: Dynamic Stiffness Matrix

[Kg] = Stiffness_Mass_matrix( Nnode,Nelement,Ielement,Coordinate,density,mu,cf);

%% Lamb dispersion
w_sca = cf/CT*h; % 无量纲频率
Fun = @lamb;
[root,~,modes] = get_wavenumber(w_sca,Fun); % 得到无量纲波数和模态个数
root = root(end:-1:1); % 按低阶到高阶排列
[Amp] = get_amplitude( root,w_sca,h,Fun );


%% Traction_Known_R (-Kg*U_in+t_401_in+t_501_in)
[ F,U_in ] = Traction_Known_R( Nnode,mode_in,f,CT,CT,h,Coordinate,Ielement,root,Kg,W,dW,dL,L_nod_s,L_nod_d,Amp,bnd_401_e,bnd_501_e,mu,mu );
% [ F,U_in ] = Traction_Known_R( Nnode,mode_in,f,CT1,he,Coordinate,Ielement,root,Kg,W,dW,dL,L_nod_s,L_nod_d,Amp,bnd_401_e,bnd_501_e,mu1 );

%% UnKnown Traction (left hand side)
[ Kg,zint_inv ] = Traction_scatter_401_501( bnd_401,bnd_501,bnd_401_e,bnd_501_e,root,f,CT,CT2,h,mu,mu,Coordinate,Ielement,dW,Amp,Kg );

U_sc = Kg\F;

for i = 1:size(bnd_401,2)
	indx_n = bnd_401(1,i);
	Usc_401(i,1)= U_sc(indx_n,1);
end

for i = 1:size(bnd_501,2)
	indx_n = bnd_501(1,i);
	Usc_501(i,1)= U_sc(indx_n,1);
end

%% Re & Tr coefficients
n_mode_max = size(root,2);
[ Re,Tr,eng_R,eng_T ] = Co_Eng_RaT( n_mode_max,mode_in,bnd_401_e,bnd_501_e,Coordinate,Ielement,root,f,CT,CT,mu,mu,h,dW,Amp,zint_inv,U_sc );
