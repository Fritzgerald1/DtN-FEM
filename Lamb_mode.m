function [ z_um,zderiv_um ] = Lamb_mode(m,Amp,i_layer,kth1,kth2,kh,x1,x2,i_pm,i_conj)
%SH_MODE 此处显示有关此函数的摘要
%   i_pm = 1(positive direction), -1(negative direction)
he= 5*10^(-4);
kth = [kth1 kth2];
zsqrtk = sqrt(kth(1,i_layer)^2-kh^2);
j = 2*(i_layer-1)+1;
z_um = (Amp(j,m)*exp(1i*zsqrtk*x2)+Amp(j+1,m)*exp(-1i*zsqrtk*x2))*...
        exp(i_pm*1i*kh*x1);
zderiv_um(1,1)=i_pm*1i*kh*z_um;  %%
zderiv_um(2,1)=1i*zsqrtk*(Amp(j,m)*exp(1i*zsqrtk*x2)-Amp(j+1,m)*exp(-1i*zsqrtk*x2))*...
                exp(i_pm*1i*kh*x1);
if i_conj == 1
    z_um = conj(z_um);
    zderiv_um(:,1) = conj(zderiv_um(:,1));
end


end

