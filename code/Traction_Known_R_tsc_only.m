function [ F ] = Traction_Known_R_tsc_only( Nnode,Nelement,Nelement_crack,L_ele,f,Coordinate,Ielement,bnd_101_e,CT1,CT2,he,root,Amp,mode_in,mu1,mu2,dW,dL )
% 101 boundary
F=spalloc(Nnode,1,10000);
Ne_101_top = Nelement-L_ele+1;
Ne_101_bot = L_ele;
Ne_101_mid_top_l = round(Nelement/2+L_ele/2-Nelement_crack/2+1);
Ne_101_mid_top_r = round(Nelement/2+L_ele/2+Nelement_crack/2);
Ne_101_mid_bot_l = round(Nelement/2-L_ele/2-Nelement_crack/2+1);
Ne_101_mid_bot_r = round(Nelement/2-L_ele/2+Nelement_crack/2);
for j = 1:size(bnd_101_e,2)
    indx_e = bnd_101_e(1,j);
    if indx_e <=Ne_101_bot || (indx_e >= Ne_101_mid_top_l && indx_e <= Ne_101_mid_top_r)
       kh = root(1,mode_in);
       kth1 = 2*pi*f/CT1*he;
       kth2 = 2*pi*f/CT2*he;
       if indx_e <= Ne_101_bot
           mu = 1;
           i_layer = 1;
       else
           mu = mu2/mu1;
           i_layer = 2;
       end
       xs(1,1) = Coordinate(Ielement(indx_e,1),1);
       xs(1,2) = Coordinate(Ielement(indx_e,1),2);
       xs(2,1) = Coordinate(Ielement(indx_e,5),1);
       xs(2,2) = Coordinate(Ielement(indx_e,5),2);
       xs(3,1) = Coordinate(Ielement(indx_e,2),1);
       xs(3,2) = Coordinate(Ielement(indx_e,2),2);
       any(1,1)= (Coordinate(Ielement(indx_e,2),2)-Coordinate(Ielement(indx_e,1),2))/dW*he;
       any(1,2)= (Coordinate(Ielement(indx_e,1),1)-Coordinate(Ielement(indx_e,2),1))/dL*he;
       [int_NT_tin ] = inwv(mode_in,mu,i_layer,Amp,kth1,kth2,kh,xs,any,101);
       F(Ielement(indx_e,1),1) = F(Ielement(indx_e,1),1) - int_NT_tin(1,1);
       F(Ielement(indx_e,5),1) = F(Ielement(indx_e,5),1) - int_NT_tin(2,1);
       F(Ielement(indx_e,2),1) = F(Ielement(indx_e,2),1) - int_NT_tin(3,1);
    elseif indx_e >=Ne_101_top || (indx_e >= Ne_101_mid_bot_l && indx_e <= Ne_101_mid_bot_r)
       kh = root(1,mode_in);
       kth1 = 2*pi*f/CT1*he;
       kth2 = 2*pi*f/CT2*he;
       if indx_e >=Ne_101_top
           mu = mu2/mu1;
           i_layer = 2;
       else
           mu = 1;
           i_layer = 1;
       end
       xs(1,1) = Coordinate(Ielement(indx_e,3),1);
       xs(1,2) = Coordinate(Ielement(indx_e,3),2);
       xs(2,1) = Coordinate(Ielement(indx_e,7),1);
       xs(2,2) = Coordinate(Ielement(indx_e,7),2);
       xs(3,1) = Coordinate(Ielement(indx_e,4),1);
       xs(3,2) = Coordinate(Ielement(indx_e,4),2);
       any(1,1)= (Coordinate(Ielement(indx_e,4),2)-Coordinate(Ielement(indx_e,3),2))/dW*he;%%
       any(1,2)= (Coordinate(Ielement(indx_e,3),1)-Coordinate(Ielement(indx_e,4),1))/dL*he;%%
       [int_NT_tin ] = inwv(mode_in,mu,i_layer,Amp,kth1,kth2,kh,xs,any,102);
       F(Ielement(indx_e,3),1) = F(Ielement(indx_e,3),1) - int_NT_tin(1,1);
       F(Ielement(indx_e,7),1) = F(Ielement(indx_e,7),1) - int_NT_tin(2,1);
       F(Ielement(indx_e,4),1) = F(Ielement(indx_e,4),1) - int_NT_tin(3,1);
   end
end

end

