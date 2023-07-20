function [ Re,Tr,eng_R,eng_T,zbeta_R,zbeta_T ] = Co_Eng_RaT( n_mode_max,mode_in,bnd_401_e,bnd_501_e,Coordinate,Ielement,root,f,CT1,CT2,mu1,mu2,he,dW,Amp,zint_inv,U_sc )
zbeta_R = zeros(n_mode_max,1);
for i1 = 1:size(bnd_401_e,2)
    indx_e = bnd_401_e(1,i1);
    xf(1,1) = (Coordinate(Ielement(indx_e,1),1)+Coordinate(Ielement(indx_e,4),1))*0.5;
    for imode = 1:n_mode_max
         for jmode = 1:n_mode_max
             xs(1,1)=0;
             xs(1,2)=Coordinate(Ielement(indx_e,1),2);
             xs(2,1)=0;
             xs(2,2)=Coordinate(Ielement(indx_e,8),2);
             xs(3,1)=0;
             xs(3,2)=Coordinate(Ielement(indx_e,4),2);
             kth1 = 2*pi*f/CT1*he;
             kth2 = 2*pi*f/CT2*he;
             kh = root(1,jmode);
             if i1 <= (size(bnd_401_e,2)/2)
                i_layer = 1;
             else
                i_layer = 2;
             end
             [ zint_fm ] = int_base_um_conjg( xs,jmode,i_layer,Amp,kth1,kth2,kh,0 );
             zbeta_R(imode,1) = zbeta_R(imode,1)+...
                                exp(1i*root(1,imode)*xf(1,1))*zint_inv(imode,jmode)*... 
                                zint_fm*[U_sc(Ielement(indx_e,1),1);U_sc(Ielement(indx_e,8),1);U_sc(Ielement(indx_e,4),1)];
        end
    end
end

zbeta_T = zeros(n_mode_max,1);
for i1 = 1:size(bnd_501_e,2)
    indx_e = bnd_501_e(1,i1);
    xf(1,1) = (Coordinate(Ielement(indx_e,2),1)+Coordinate(Ielement(indx_e,3),1))*0.5;
    for imode = 1:n_mode_max
         for jmode = 1:n_mode_max
             xs(1,1)=0;
             xs(1,2)=Coordinate(Ielement(indx_e,2),2);
             xs(2,1)=0;
             xs(2,2)=Coordinate(Ielement(indx_e,6),2);
             xs(3,1)=0;
             xs(3,2)=Coordinate(Ielement(indx_e,3),2);
             kth1 = 2*pi*f/CT1*he;
             kth2 = 2*pi*f/CT2*he;
             kh = root(1,jmode);
             if i1 <= (size(bnd_501_e,2)/2)
                i_layer = 1;
             else
                i_layer = 2;
             end
             [ zint_fm ] = int_base_um_conjg( xs,jmode,i_layer,Amp,kth1,kth2,kh,1 );
             zbeta_T(imode,1) = zbeta_T(imode,1)+...
                                exp(-1i*root(1,imode)*xf(1,1))*zint_inv(imode,jmode)*...
                                zint_fm*[U_sc(Ielement(indx_e,2),1);U_sc(Ielement(indx_e,6),1);U_sc(Ielement(indx_e,3),1)];
           
        end
    end
end

for imode=1:n_mode_max
    kth1 = 2*pi*f/CT1*he;
    kth2 = 2*pi*f/CT2*he;
    kih = root(1,imode);
    n_layer = 2;
    miu = [mu1 mu2];
    [ Pm ] = power_eval_t_conju( imode,imode,n_layer,miu,kth1,kth2,kih,kih,Amp,0,he,1 );
    pm(imode,1)=Pm;
end


if mode_in <= n_mode_max
    eng_R = 0;
    eng_T = 0;
    for imode=1:n_mode_max
        eng_R =eng_R+pm(imode,1)/pm(mode_in,1)*(conj(zbeta_R(imode,1))*zbeta_R(imode));
        if imode == mode_in
            Re(imode,1)=abs(zbeta_R(imode,1));
            Tr(imode,1)=abs(1+zbeta_T(imode,1));
            eng_T=eng_T+pm(imode,1)/pm(mode_in,1)*(conj(1+zbeta_T(imode,1))*(1+zbeta_T(imode,1)));
        else
            Re(imode,1)=abs(zbeta_R(imode,1));
            Tr(imode,1)=abs(zbeta_T(imode,1));
            eng_T=eng_T+pm(imode,1)/pm(mode_in,1)*(conj(zbeta_T(imode,1))*(zbeta_T(imode,1)));
        end
    end
end


end

