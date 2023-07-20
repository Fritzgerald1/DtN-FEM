function [ zint_u_mode ] = int_u_dz(imode,jmode,n_layer,kth1,kth2,kih,kjh,zamp_mode,x1,h)
% double layer
y = [-1 0 1];
kth = [kth1 kth2];
zint_u_mode = 0;
%for i=1:n_layer
%    y2=y(1,i+1);
%    y1=y(1,i);
%    zsqrtki = sqrt(kth(1,i)^2-kih^2);
%    zsqrtkj = sqrt(kth(1,i)^2-kjh^2);
%    for j=1:2
%        for k=1:2
%            zdum1 = conj((-1)^(j-1)*1i*zsqrtki) + (-1)^(k-1)*1i*zsqrtkj;
%            if abs(zdum1) < 1*10^(-8)
%                zdum2 = y2-y1+0.5*zdum1*(y2^2-y1^2)+1/6*zdum1*zdum1*(y2^3-y1^3);
%            else
%                zdum2 = (exp(zdum1*y2)-exp(zdum1*y1))/zdum1;
%            end
%            zint_u_mode = zint_u_mode + conj(zamp_mode(j+2*(i-1),imode))*zamp_mode(k+2*(i-1),jmode)*zdum2*...
%                          exp(1i*(-kih+kjh)*x1);
%        end
%    end
%end

n = 8;
[x,w]=gauss_integration(n);
for i =1:n_layer
    y2=y(1,i+1);
    y1=y(1,i);
    qq = 0.5*(y2-y1);
    zsqrtki = sqrt(kth(1,i)^2-kih^2);
    zsqrtkj = sqrt(kth(1,i)^2-kjh^2);
    for iw=1:n
        ys1=x1;
        ys2=((y2-y1)*x(iw)+y1+y2)*0.5;
        z_um1=(zamp_mode(2*i-1,imode)*exp(1i*zsqrtki*ys2)+zamp_mode(2*i,imode)*exp(-1i*zsqrtki*ys2))*exp(1i*kih*ys1);
        z_um2=(zamp_mode(2*i-1,jmode)*exp(1i*zsqrtkj*ys2)+zamp_mode(2*i,jmode)*exp(-1i*zsqrtkj*ys2))*exp(1i*kjh*ys1);
        z_um1=conj(z_um1);
        zint_u_mode = zint_u_mode+z_um1*z_um2*qq*w(iw);
    end
end

end