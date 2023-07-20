function [ Kg,zint_inv,zint_amp_mode] = Traction_scatter_401_501( bnd_401,bnd_501,bnd_401_e,bnd_501_e,root,f,CT1,CT2,he,mu1,mu2,Coordinate,Ielement,dW,Amp,Kg )
% 401 and 501 boundary (scatter component) %
% int_NT_tm
L_int_NTtm = zeros(size(bnd_401,2),size(root,2));
R_int_NTtm = zeros(size(bnd_501,2),size(root,2));
for i = 1:size(root,2)
   for i1 = 1:size(bnd_401_e,2)
       indx_e = bnd_401_e(1,i1);
       kh = root(1,i);
       kth1 = 2*pi*f/CT1*he;
       kth2 = 2*pi*f/CT2*he;
       if i1 <= (size(bnd_401_e,2)/2)
            i_layer = 1;
            mu = 1;
       else
            i_layer = 2;
            mu = mu2/mu1;
       end
       xs(1,1) = 0;
       xs(1,2) = Coordinate(Ielement(indx_e,1),2);
       xs(2,1) = 0;
       xs(2,2) = Coordinate(Ielement(indx_e,8),2);
       xs(3,1) = 0;
       xs(3,2) = Coordinate(Ielement(indx_e,4),2);
       [ zint_NTtm ] = int_NT_tm( i,mu,i_layer,xs,Amp,kth1,kth2,kh,0 );
       L_int_NTtm(2*i1-1,i) = L_int_NTtm(2*i1-1,i)+zint_NTtm(1,1);  %%
       L_int_NTtm(2*i1,i) = L_int_NTtm(2*i1,i)+zint_NTtm(2,1); %%
       L_int_NTtm(2*i1+1,i) = L_int_NTtm(2*i1+1,i)+zint_NTtm(3,1); %%
   end
   for i1 = 1:size(bnd_501_e,2)
       indx_e = bnd_501_e(1,i1);
       kh = root(1,i);
       kth1 = 2*pi*f/CT1*he;
       kth2 = 2*pi*f/CT2*he;
       if i1 <= (size(bnd_501_e,2)/2)
            i_layer = 1;
            mu = 1;
       else
            i_layer = 2;
            mu = mu2/mu1;
       end
       xs(1,1) = 0;
       xs(1,2) = Coordinate(Ielement(indx_e,2),2);
       xs(2,1) = 0;
       xs(2,2) = Coordinate(Ielement(indx_e,6),2);
       xs(3,1) = 0;
       xs(3,2) = Coordinate(Ielement(indx_e,3),2);
       [ zint_NTtm ] = int_NT_tm( i,mu,i_layer,xs,Amp,kth1,kth2,kh,1 );
       R_int_NTtm(2*i1-1,i) = R_int_NTtm(2*i1-1,i)+zint_NTtm(1,1); %%
       R_int_NTtm(2*i1,i) = R_int_NTtm(2*i1,i)+zint_NTtm(2,1); %%
       R_int_NTtm(2*i1+1,i) = R_int_NTtm(2*i1+1,i)+zint_NTtm(3,1); %%
   end
end
   
% Dmn^(-1)
n_mode_max = size(root,2);
for i = 1:n_mode_max
    for j = 1:n_mode_max
        kth1 = 2*pi*f/CT1*he;
        kth2 = 2*pi*f/CT2*he;
        n_layer = 2;
        kih = root(1,i);
        kjh = root(1,j);
        zint_u_mode = int_u_dz(i,j,n_layer,kth1,kth2,kih,kjh,Amp,0,he);%%%
        zint_amp_mode(i,j) = zint_u_mode;
    end
end

ztmp_vec = zeros(n_mode_max,1);
for j=1:n_mode_max
    ztmp_vec(:,1) = 0;
    ztmp_vec(j,1) = 1;
    x = zint_amp_mode\ztmp_vec;
    for i=1:n_mode_max
        zint_inv(i,j) = x(i,1);
    end
end

% int_base_um_conj
%jindx_leftedge = 0;
%jindx_rightedge = 0;
zintu_leftedge = zeros(n_mode_max,size(bnd_401,2));
zintu_rightedge = zeros(n_mode_max,size(bnd_501,2));
for imode = 1:n_mode_max
    for j = 1:size(bnd_401_e,2)
        indx_e = bnd_401_e(1,j);   
        xs(1,1)=0;
        xs(1,2)=Coordinate(Ielement(indx_e,1),2);
        xs(2,1)=0;
        xs(2,2)=Coordinate(Ielement(indx_e,8),2);
        xs(3,1)=0;
        xs(3,2)=Coordinate(Ielement(indx_e,4),2);
        kth1 = 2*pi*f/CT1*he;
        kth2 = 2*pi*f/CT2*he;
        kh = root(1,imode);
        if j <= (size(bnd_401_e,2)/2)
            i_layer = 1;
        else
            i_layer = 2;
        end
        [ zint_fm ] = int_base_um_conjg( xs,imode,i_layer,Amp,kth1,kth2,kh,0 );
        zintu_leftedge(imode,2*j-1) = zintu_leftedge(imode,2*j-1)+zint_fm(1,1);
        zintu_leftedge(imode,2*j) = zintu_leftedge(imode,2*j)+zint_fm(1,2);
        zintu_leftedge(imode,2*j+1) = zintu_leftedge(imode,2*j+1)+zint_fm(1,3);
    end
end

for imode = 1:n_mode_max
    for j = 1:size(bnd_501_e,2)
        indx_e = bnd_501_e(1,j);
        xs(1,1)=0;
        xs(1,2)=Coordinate(Ielement(indx_e,2),2);
        xs(2,1)=0;
        xs(2,2)=Coordinate(Ielement(indx_e,6),2);
        xs(3,1)=0;
        xs(3,2)=Coordinate(Ielement(indx_e,3),2);
        kth1 = 2*pi*f/CT1*he;
        kth2 = 2*pi*f/CT2*he;
        kh = root(1,imode);
        if j <= (size(bnd_501_e,2)/2)
            i_layer = 1;
        else
            i_layer = 2;
        end
        [ zint_fm ] = int_base_um_conjg( xs,imode,i_layer,Amp,kth1,kth2,kh,1 );
        zintu_rightedge(imode,2*j-1) = zintu_rightedge(imode,2*j-1)+zint_fm(1,1);
        zintu_rightedge(imode,2*j) = zintu_rightedge(imode,2*j)+zint_fm(1,2);
        zintu_rightedge(imode,2*j+1) = zintu_rightedge(imode,2*j+1)+zint_fm(1,3);
    end
end        

% combine int_NTtm*zint_inv*zintu_leftedge/zintu_rightedge
L_int = 0;
R_int = 0;
for imode = 1:n_mode_max
    for jmode = 1:n_mode_max
        L_int = L_int+L_int_NTtm(:,imode)*zint_inv(imode,jmode)*zintu_leftedge(jmode,:);
        R_int = R_int+R_int_NTtm(:,imode)*zint_inv(imode,jmode)*zintu_rightedge(jmode,:);
    end
end

% 401 boundary
for i = 1:size(bnd_401,2)
    indx_n_i = bnd_401(1,i);
    for j = 1:size(bnd_401,2)
        indx_n_j = bnd_401(1,j);
        Kg(indx_n_i,indx_n_j)=Kg(indx_n_i,indx_n_j)-L_int(i,j); %%
    end
end
    
% 501 boundary
for i = 1:size(bnd_501,2)
    indx_n_i = bnd_501(1,i);
    for j = 1:size(bnd_501,2)
        indx_n_j = bnd_501(1,j);
        Kg(indx_n_i,indx_n_j)=Kg(indx_n_i,indx_n_j)-R_int(i,j); %%
    end
end

end

