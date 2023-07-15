function [ F,U_in ] = Traction_Known_R( Nnode,mode_in,f,CT1,CT2,he,Coordinate,Ielement,root,Kg,W,dW,dL,L_nod_s,L_nod_d,Amp,bnd_401_e,bnd_501_e,mu1,mu2 )
% 101,401,501 boundary  (incident component Kg*U_in)
U_in=zeros(Nnode,1);
F=zeros(Nnode,1);
for j = 1:Nnode
    kh = root(1,mode_in);
    kth1 = 2*pi*f/CT1*he;
    kth2 = 2*pi*f/CT2*he;
    x = Coordinate(j,1);
    y = Coordinate(j,2);
    if j <= ((fix(W/dW)/2+1)*L_nod_s+fix(W/dW)/2*L_nod_d)
        i_layer = 1;
    else
        i_layer = 2;
    end
    [ z_um,zderiv_um] = Lamb_mode(mode_in,Amp,i_layer,kth1,kth2,kh,x,y,1,0);
    U_in(j,1) = z_um;
end

F = F-Kg*U_in;

% 401,501 boundary  (incident component int_NT_tin)
for j = 1:size(bnd_401_e,2)
    indx_e = bnd_401_e(1,j);
    kh = root(1,mode_in);
    kth1 = 2*pi*f/CT1*he;
    kth2 = 2*pi*f/CT2*he;
    if j <= (size(bnd_401_e,2)/2)
        i_layer = 1;
        mu = 1;
    else
        i_layer = 2;
        mu = mu2/mu1;
	end
    xs(1,1) = Coordinate(Ielement(indx_e,4),1);
    xs(1,2) = Coordinate(Ielement(indx_e,4),2);
    xs(2,1) = Coordinate(Ielement(indx_e,8),1);
    xs(2,2) = Coordinate(Ielement(indx_e,8),2);
    xs(3,1) = Coordinate(Ielement(indx_e,1),1);
    xs(3,2) = Coordinate(Ielement(indx_e,1),2);
    any(1,1)= (Coordinate(Ielement(indx_e,1),2)-Coordinate(Ielement(indx_e,4),2))/dW*he; 
    any(1,2)= (Coordinate(Ielement(indx_e,4),1)-Coordinate(Ielement(indx_e,1),1))/dL*he;
    [int_NT_tin ] = inwv(mode_in,mu,i_layer,Amp,kth1,kth2,kh,xs,any,401);
    F(Ielement(indx_e,4),1) = F(Ielement(indx_e,4),1) + int_NT_tin(1,1); %%
    F(Ielement(indx_e,8),1) = F(Ielement(indx_e,8),1) + int_NT_tin(2,1); %%
    F(Ielement(indx_e,1),1) = F(Ielement(indx_e,1),1) + int_NT_tin(3,1); 
end

for j = 1:size(bnd_501_e,2)
    indx_e = bnd_501_e(1,j);
    kh = root(1,mode_in);
    kth1 = 2*pi*f/CT1*he;
    kth2 = 2*pi*f/CT2*he;
    if j <= (size(bnd_501_e,2)/2)
        i_layer = 1;
        mu = 1;
    else
        i_layer = 2;
        mu = mu2/mu1;
    end
    xs(1,1) = Coordinate(Ielement(indx_e,2),1);
    xs(1,2) = Coordinate(Ielement(indx_e,2),2);
    xs(2,1) = Coordinate(Ielement(indx_e,6),1);
    xs(2,2) = Coordinate(Ielement(indx_e,6),2);
    xs(3,1) = Coordinate(Ielement(indx_e,3),1);
    xs(3,2) = Coordinate(Ielement(indx_e,3),2);
    any(1,1)= (Coordinate(Ielement(indx_e,3),2)-Coordinate(Ielement(indx_e,2),2))/dW*he;
    any(1,2)= (Coordinate(Ielement(indx_e,2),1)-Coordinate(Ielement(indx_e,3),1))/dL*he;
    [int_NT_tin ] = inwv(mode_in,mu,i_layer,Amp,kth1,kth2,kh,xs,any,501);
    F(Ielement(indx_e,2),1) = F(Ielement(indx_e,2),1) + int_NT_tin(1,1);  %%
    F(Ielement(indx_e,6),1) = F(Ielement(indx_e,6),1) + int_NT_tin(2,1); %%
    F(Ielement(indx_e,3),1) = F(Ielement(indx_e,3),1) + int_NT_tin(3,1);
end


end

