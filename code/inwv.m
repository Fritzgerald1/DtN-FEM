function [int_NT_tin ] = inwv(m,mu,i_layer,Amp,kth1,kth2,kh,xs,any,i_symbol)
int_NT_tin = zeros(3,1);
rs = sqrt((xs(3,1)-xs(1,1))^2+(xs(3,2)-xs(1,2))^2);
qq = 0.5*rs;
n = 8;
he = 5*10^(-4);
[x,w]=gauss_integration(n);
for iw = 1:n
    ys1 =((xs(3,1)-xs(1,1))*x(iw)+xs(3,1)+xs(1,1))*0.5;
    ys2 =((xs(3,2)-xs(1,2))*x(iw)+xs(3,2)+xs(1,2))*0.5;
    [ z_um,zderiv_um ] = SH_mode(m,Amp,i_layer,kth1,kth2,kh,ys1,ys2,1,0);
    if i_symbol == 401   % 4-8-1
        NT(1,1) = 1/(rs*qq)*(xs(3,2)-ys2)*(xs(2,2)-ys2);
        NT(2,1) = -1/(qq^2)*(xs(3,2)-ys2)*(xs(1,2)-ys2);
        NT(3,1) = 1/(rs*qq)*(xs(2,2)-ys2)*(xs(1,2)-ys2);
    elseif i_symbol == 501  % 2-6-3
        NT(1,1) = 1/(rs*qq)*(xs(2,2)-ys2)*(xs(3,2)-ys2);
        NT(2,1) = -1/(qq^2)*(xs(1,2)-ys2)*(xs(3,2)-ys2);
        NT(3,1) = 1/(rs*qq)*(xs(1,2)-ys2)*(xs(2,2)-ys2);
    elseif i_symbol == 101  % 1-5-2
        NT(1,1) = 1/(rs*qq)*(xs(2,1)-ys1)*(xs(3,1)-ys1);
        NT(2,1) = -1/(qq^2)*(xs(1,1)-ys1)*(xs(3,1)-ys1);
        NT(3,1) = 1/(rs*qq)*(xs(1,1)-ys1)*(xs(2,1)-ys1);
    elseif i_symbol == 102  % 3-7-4
        NT(1,1) = 1/(rs*qq)*(xs(2,1)-ys1)*(xs(3,1)-ys1);
        NT(2,1) = -1/(qq^2)*(xs(1,1)-ys1)*(xs(3,1)-ys1); 
        NT(3,1) = 1/(rs*qq)*(xs(1,1)-ys1)*(xs(2,1)-ys1); 
    end
    tin = (any(1,1)*zderiv_um(1,1)+any(1,2)*zderiv_um(2,1));
    int_NT_tin = int_NT_tin+NT*mu*tin*qq*w(iw);
end
end

