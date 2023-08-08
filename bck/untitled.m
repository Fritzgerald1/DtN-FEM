density = 2.7e3; 
mu = 26e9; 

f = 1.5*10^6; % 入射频率
cf = 2*pi*f;
%% 
% 

[Nnode, Nelement, Coordinate, Ielement, dL, dW] = read_mesh_info();

[K,M] = Stiffness_Mass_matrix( Nnode,Nelement,Ielement,Coordinate,density,mu,cf);