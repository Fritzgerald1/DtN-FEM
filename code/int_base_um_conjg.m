function [ zint_fm ] = int_base_um_conjg( xs,m,i_layer,Amp,kth1,kth2,kh,i_symbol )
zint_fm = zeros(1,3);
rs = sqrt((xs(3,1)-xs(1,1))^2+(xs(3,2)-xs(1,2))^2);
qq = 0.5*rs;
n = 8;
[x,w]=gauss_integration(n);
for iw = 1:n
    ys1 =((xs(3,1)-xs(1,1))*x(iw)+xs(3,1)+xs(1,1))*0.5;
    ys2 =((xs(3,2)-xs(1,2))*x(iw)+xs(3,2)+xs(1,2))*0.5;
    [ z_um,zderiv_um ] = SH_mode(m,Amp,i_layer,kth1,kth2,kh,ys1,ys2,1,1);
    if i_symbol == 0  %%1-8-4
        N(1,1) = 1/(qq*rs)*(xs(2,2)-ys2)*(xs(3,2)-ys2);
        N(1,2) = -1/(qq^2)*(xs(1,2)-ys2)*(xs(3,2)-ys2);
        N(1,3) = 1/(qq*rs)*(xs(1,2)-ys2)*(xs(2,2)-ys2);
    elseif i_symbol == 1
        N(1,1) = 1/(qq*rs)*(xs(2,2)-ys2)*(xs(3,2)-ys2);
        N(1,2) = -1/(qq^2)*(xs(1,2)-ys2)*(xs(3,2)-ys2);
        N(1,3) = 1/(qq*rs)*(xs(1,2)-ys2)*(xs(2,2)-ys2);
    end    
    zint_fm = zint_fm+z_um*N*qq*w(iw);
end


end

