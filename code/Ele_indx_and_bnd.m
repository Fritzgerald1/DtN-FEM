function [ bnd_401_e,bnd_501_e,bnd_101_e,bnd_401,bnd_501,bnd_101,Ielement ] = Ele_indx_and_bnd( L_ele,L_nod_s,L_nod_d,Nelement,Nelement_crack,Nnode,Nnode_crack)
icount=1;
cNode=1;
count_401_e = 0;
count_501_e = 0;
count_101_e = 0;
count_401 = 1;
count_501 = 1;
count_101 = 1;
Ne_mid_f = 0;
D_Lcount = 0;
Ne_101_top = Nelement-L_ele+1;
Ne_101_bot = L_ele;
Ne_101_mid_top_l = round(Nelement/2+L_ele/2-Nelement_crack/2+1);
Ne_101_mid_top_r = round(Nelement/2+L_ele/2+Nelement_crack/2);
Ne_101_mid_bot_l = round(Nelement/2-L_ele/2-Nelement_crack/2+1);
Ne_101_mid_bot_r = round(Nelement/2-L_ele/2+Nelement_crack/2);
L_nod = L_nod_s+L_nod_d;
Ielement=zeros(Nelement,8);
for Ne=1:Nelement
    if Ne == Ne_101_mid_top_l
        Ielement(Ne,:)=[cNode (Nnode-Nnode_crack+2) cNode+L_nod+2 cNode+L_nod (Nnode-Nnode_crack+1) cNode+L_nod_s-D_Lcount+1 cNode+L_nod+1 cNode+L_nod_s-D_Lcount];
        Ne_mid_f = Ne;
    elseif Ne == Ne_101_mid_top_r
        Ielement(Ne,:)=[Nnode-1 cNode+2 cNode+L_nod+2 cNode+L_nod Nnode cNode+L_nod_s-D_Lcount+1 cNode+L_nod+1 cNode+L_nod_s-D_Lcount];
    elseif  Ne > Ne_101_mid_top_l && Ne < Ne_101_mid_top_r
        Ielement(Ne,:)=[(Nnode-Nnode_crack+2*(Ne-Ne_mid_f)) (Nnode-Nnode_crack+2*(Ne-Ne_mid_f+1)) cNode+L_nod+2 cNode+L_nod (Nnode-Nnode_crack+2*(Ne-Ne_mid_f)+1) cNode+L_nod_s+1-D_Lcount cNode+L_nod+1 cNode+L_nod_s-D_Lcount];
    else
        Ielement(Ne,:)=[cNode cNode+2 cNode+L_nod+2 cNode+L_nod cNode+1 cNode+L_nod_s-D_Lcount+1 cNode+L_nod+1 cNode+L_nod_s-D_Lcount];
    end
    
    cNode=cNode+2;
    D_Lcount = D_Lcount +1;
    if mod(Ne,L_ele)==0
        icount=icount+1;
        cNode=(icount-1)*L_nod+1;
        D_Lcount = 0;
    end
    %%% bnd_401, bnd_501, bnd_101
    if mod(Ne,L_ele) == 1
        count_401_e = count_401_e+1;
        bnd_401_e(1,count_401_e) = Ne;
        if count_401 == 1
            bnd_401(1,count_401) = Ielement(Ne,1);
            bnd_401(1,count_401+1) = Ielement(Ne,8);
            bnd_401(1,count_401+2) = Ielement(Ne,4);
        else
            bnd_401(1,count_401+1) = Ielement(Ne,8);
            bnd_401(1,count_401+2) = Ielement(Ne,4);
        end
        count_401 = count_401+2;
    elseif mod(Ne,L_ele) == 0
        count_501_e = count_501_e+1;
        bnd_501_e(1,count_501_e) = Ne;
        if count_501 == 1
            bnd_501(1,count_501) = Ielement(Ne,2);
            bnd_501(1,count_501+1) = Ielement(Ne,6);
            bnd_501(1,count_501+2) = Ielement(Ne,3);
        else
            bnd_501(1,count_501+1) = Ielement(Ne,6);
            bnd_501(1,count_501+2) = Ielement(Ne,3);
        end
        count_501 = count_501+2;
    end
    if Ne <= Ne_101_bot || Ne >= Ne_101_top
        count_101_e = count_101_e+1;
        bnd_101_e(1,count_101_e) = Ne;
        if Ne == Ne_101_bot
            bnd_101(1,count_101) = Ielement(Ne,1);
            bnd_101(1,count_101+1) = Ielement(Ne,5);
            bnd_101(1,count_101+2) = Ielement(Ne,2);
            count_101 = count_101+3;
        elseif Ne < Ne_101_bot
            bnd_101(1,count_101) = Ielement(Ne,1);
            bnd_101(1,count_101+1) = Ielement(Ne,5);
            count_101 = count_101+2;
        elseif Ne == Ne_101_top
            bnd_101(1,count_101) = Ielement(Ne,4);
            bnd_101(1,count_101+1) = Ielement(Ne,7);
            bnd_101(1,count_101+2) = Ielement(Ne,3);
            count_101 = count_101+3;
        elseif Ne > Ne_101_top
            bnd_101(1,count_101) = Ielement(Ne,7);
            bnd_101(1,count_101+1) = Ielement(Ne,3);
            count_101 = count_101+2;
        end         
    elseif Ne >= Ne_101_mid_top_l && Ne <= Ne_101_mid_top_r
        count_101_e = count_101_e+1;
        bnd_101_e(1,count_101_e) = Ne;
        if Ne < Ne_101_mid_top_r
            bnd_101(1,count_101) = Ielement(Ne,5);
            bnd_101(1,count_101+1) = Ielement(Ne,2);
            count_101 = count_101+2;  
        elseif Ne == Ne_101_mid_top_r
            bnd_101(1,count_101) = Ielement(Ne,5);
            count_101 = count_101+1;
        end    
    elseif Ne >= Ne_101_mid_bot_l && Ne <=Ne_101_mid_bot_r
        count_101_e = count_101_e+1;
        bnd_101_e(1,count_101_e) = Ne;
        if Ne == Ne_101_mid_bot_l
            bnd_101(1,count_101) = Ielement(Ne,4);
            bnd_101(1,count_101+1) = Ielement(Ne,7);
            bnd_101(1,count_101+2) = Ielement(Ne,3);
            count_101 = count_101+3;
        elseif Ne == Ne_101_mid_bot_r
            bnd_101(1,count_101) = Ielement(Ne,7);
            bnd_101(1,count_101+1) = Ielement(Ne,3);
            count_101 = count_101+2;    
        else
            bnd_101(1,count_101) = Ielement(Ne,7);
            bnd_101(1,count_101+1) = Ielement(Ne,3);
            count_101 = count_101+2;
        end  
    end
end

end

