clc;
clear;
% input settings %
mlayer = 2;
% material parameter %
mu1 = 79*10^3;
mu2 = 79*10^3;
density1 = 7.8*10^3;
density2 = 7.8*10^3;
CT1 = sqrt(mu1/density1);
CT2 = sqrt(mu2/density2);
tol = 1*10^(-12);
df = 0.1;
% calculation %


for inf = 0:100
    f = inf * df;
    sum = 0;
    khmin = 0;
    khmax = 2*pi*f;
    dkh = (khmin-khmax)/10000;
    for ikh = 1:10000
        kth1 = 2*pi*f/CT1*0.5;
        kth2 = 2*pi*f/CT2*0.5;
        a = khmax + (ikh-2)*dkh;
        b = a + dkh;
        c = a + 2*dkh;
        [det_a] = double_fk(a,kth1,kth2,mu1,mu2);
        [det_b] = double_fk(b,kth1,kth2,mu1,mu2);
        [det_c] = double_fk(c,kth1,kth2,mu1,mu2);
        if ikh > 1 && inf > 0
            if det_b < det_a && det_b < det_c
                [xmin,error] = brent_fk_double(a,b,c,tol,kth1,kth2,mu1,mu2);
                 if error < -5
                    sum= sum +1;
                    root(inf,sum)=xmin;
                    root_d(inf,sum)=error;
                end
            end
        end   
    end 
end


for in_f = 1:size(root,1)
    for n_r = 1:size(root,2)
        if root(in_f,n_r) ~= 0
            kth1 = 2*pi*in_f*df/CT1*0.5;
            kth2 = 2*pi*in_f*df/CT2*0.5;
            [A,D] = double_fk_amplitude(root(in_f,n_r),kth1,kth2,mu1,mu2 );
            Amp((4*in_f-3):4*in_f,n_r) = A;
        else
            break
        end
    end
end

for in_f = 32
    for n_r = 3
        for ih = 0:10
            h = -0.5+0.1*ih;
            if ih <= 5
                kt = 2*pi*in_f*0.1/CT1;
                k = root(in_f,n_r)/0.5;
                alpha = sqrt((kt/k)^2-1);
                u(ih+1) = Amp((4*in_f-3),n_r)*exp(1i*k*alpha*h)+Amp((4*in_f-2),n_r)*exp(1i*k*-alpha*h);
            else
                kt = 2*pi*in_f*0.1/CT2;
                k = root(in_f,n_r)/0.5;
                alpha = sqrt((kt/k)^2-1);
                u(ih+1) = Amp((4*in_f-1),n_r)*exp(1i*k*alpha*h)+Amp((4*in_f),n_r)*exp(1i*k*-alpha*h);
            end
        end
    end
end
                
h = -0.5:0.1:0.5;
figure
plot(u,h);



%f=df:df:10;%
%for i3=1:size(root,2)
 %   for i4=1:size(f,2)
  %      if root(i4,i3)~=0
  %          ft=i4*df; %
  %          break
   %     end
  %  end
 %   fi=ft:df:10; %
%    plot(fi,root(i4:size(f,2),i3));
%     hold on
%    axis([0 10 0 20]);
%end
                
