clc;
clear;
tic
% FEM SH wave %
% Input %
% density:matterial parameter %
% h:unit length of x3 direction %
% mu:shear modulus %
% Cf:cicular frequency %
density1=7.8*10^3;
density2=2.7*10^3;
mu1=79*10^9;
mu2=26.1*10^9;
h=1;
he = 0.5*10^(-3);
fre=2*10^6:0.025*10^6:4.5*10^6;
eng_R = zeros(1,size(fre,2));
eng_T = zeros(1,size(fre,2));
Re = zeros(5,size(fre,2));
Tr = zeros(5,size(fre,2));
%for kk = 1:size(fre,2)
%f = fre(1,kk)
%Cf = 2*pi*f;
mode_in = 2;
L=1*10^(-2);
W=1*10^(-3);
dL=L/200;
dW=W/20;
% Element settings %
% L:length of Plate, W:width of Plate %
% dL:length of element, dW:width of element %
% Nelement: number of elements, Nnode:number of Nodes %
% Ielement: Index of every nodes in local element %
L_crack = W;
A = dL*dW;
Nelement=round(L/dL*W/dW);
Nelement_crack = round(L_crack/dL);
Nnode_crack = 2*Nelement_crack-1;
Nnode=round((2*L/dL+1)*(2*W/dW+1)+Nnode_crack-fix((2*W/dW+1)/2)*L/dL); % crack area added nodes

% Node settings %
% Coordinate: coordinate of nodes %
[ Coordinate ] = Node_setting( Nnode,Nnode_crack,L_crack,L,dL,dW  );

L_ele = round(L/dL);
L_nod_s = 2*L_ele+1;
L_nod_d = L_ele+1;
[ bnd_401_e,bnd_501_e,bnd_101_e,bnd_401,bnd_501,bnd_101,Ielement ] = Ele_indx_and_bnd( L_ele,L_nod_s,L_nod_d,Nelement,Nelement_crack,Nnode,Nnode_crack);
Coordinate = Coordinate/he;

for kk = 1:size(fre,2)
f = fre(1,kk)
Cf = 2*pi*f;
% Stiffness matrix %
% ngq:integration points, total integration points is 8*8 %
% xgq:si or ti, wgq:wi or wj %
% K: global stiffness matrix %
% M: global mass matrix %
% Ke: local stiffness matrix, Me: local Mass matrix % 
[ Kg,K,M ] = Stiffness_Mass_matrix( Nnode,Nelement,Ielement,Coordinate,mu1,mu2,h,density1,density2,he,Cf );
Kg_o = Kg;
% SH dispersion 
CT1 = sqrt(mu1/density1);
CT2 = sqrt(mu2/density2);
[ root,root_d,zamp_mode,Amp ] = SH_dispersion( mu1,mu2,density1,density2,f,he );

% Traction_Known_R (-Kg*U_in+t_401_in+t_501_in)
[ F,U_in ] = Traction_Known_R( Nnode,mode_in,f,CT1,CT2,he,Coordinate,Ielement,root,Kg,W,dW,dL,L_nod_s,L_nod_d,Amp,bnd_401_e,bnd_501_e,mu1,mu2 );

% UnKnown Traction (left hand side)
[ Kg,zint_inv ] = Traction_scatter_401_501( bnd_401,bnd_501,bnd_401_e,bnd_501_e,root,f,CT1,CT2,he,mu1,mu2,Coordinate,Ielement,dW,Amp,Kg );

U_sc = Kg\F;

for i = 1:size(bnd_401,2)
    indx_n = bnd_401(1,i);
    Usc_401(i,1)= U_sc(indx_n,1);
end

for i = 1:size(bnd_501,2)
    indx_n = bnd_501(1,i);
    Usc_501(i,1)= U_sc(indx_n,1);
end

% Re,Tr,eng_R,eng_T
n_mode_max = size(root,2);
[ R,T,e_R,e_T ] = Co_Eng_RaT( n_mode_max,mode_in,bnd_401_e,bnd_501_e,Coordinate,Ielement,root,f,CT1,CT2,mu1,mu2,he,dW,Amp,zint_inv,U_sc);
for i = 1:size(R,1)
    Re(i,kk) = R(i,1);
end

for j = 1:size(T,1)
    Tr(j,kk) = T(j,1);
end
eng_R(1,kk) = e_R;
eng_T(1,kk) = e_T;
eng = e_R+e_T
end
toc