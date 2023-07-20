clc;
clear;
tic
% FEM SH wave %
% Input %
% density:matterial parameter %
% h:unit length of x3 direction %
% mu:shear modulus %
% Cf:cicular frequency %
density1=2.7*10^3;
density2=7.8*10^3;
mu1=26.1*10^9;
mu2=79*10^9;
h=1;
he = 0.5*10^(-3);
f=2*10^6;
Cf = 2*pi*f;
mode_in = 2;
L=1.5*10^(-2);
W=1*10^(-3);
dL=L/600;
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