function [ zint_NTtm ] = int_NT_tm( m,mu,i_layer,xs,Amp,kth1,kth2,kh,i_symbol )
%INT_NT_TM 此处显示有关此函数的摘要
%   此处显示详细说明
zint_NTtm = zeros(3,1);
rs = sqrt((xs(3,1)-xs(1,1))^2+(xs(3,2)-xs(1,2))^2);
qq = 0.5*rs;
n = 8;
he = 5*10^(-4);
[x,w]=gauss_integration(n);
for iw = 1:n
    ys1 =((xs(3,1)-xs(1,1))*x(iw)+xs(3,1)+xs(1,1))*0.5;
    ys2 =((xs(3,2)-xs(1,2))*x(iw)+xs(3,2)+xs(1,2))*0.5;
    %ys2 =(xs(2,2)+xs(1,2))*0.5;
    [ z_um,zderiv_um ] = SH_mode(m,Amp,i_layer,kth1,kth2,kh,ys1,ys2,1,0);
    if i_symbol == 0  %% 1-8-4
        NT(1,1) = 1/(qq*rs)*(xs(2,2)-ys2)*(xs(3,2)-ys2);
        NT(2,1) = -1/(qq^2)*(xs(1,2)-ys2)*(xs(3,2)-ys2);
        NT(3,1) = 1/(qq*rs)*(xs(1,2)-ys2)*(xs(2,2)-ys2);
    elseif i_symbol == 1 %% 2-6-3
        NT(1,1) = 1/(qq*rs)*(xs(2,2)-ys2)*(xs(3,2)-ys2);
        NT(2,1) = -1/(qq^2)*(xs(1,2)-ys2)*(xs(3,2)-ys2);
        NT(3,1) = 1/(qq*rs)*(xs(1,2)-ys2)*(xs(2,2)-ys2);
    end   
    zint_NTtm = zint_NTtm+NT*mu*1i*kh*z_um*qq*w(iw);  %%
end

end

