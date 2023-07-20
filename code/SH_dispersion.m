function [ root,root_d,zamp_mode,Amp ] = SH_dispersion( mu1,mu2,density1,density2,f,he )
CT1 = sqrt(mu1/density1);
CT2 = sqrt(mu2/density2);
tol = 1*10^(-12);
sum = 0;
cmin = 500;
khmin = 0;
khmax = 2*pi*f*he/cmin;
dkh = (khmin-khmax)/10000;
for ikh = 1:10000
    kth1 = 2*pi*f/CT1*he;
    kth2 = 2*pi*f/CT2*he;
    a = khmax + (ikh-2)*dkh;
    b = a + dkh;
    c = a + 2*dkh;
    [det_a] = double_fk(a,kth1,kth2,mu1,mu2);
    [det_b] = double_fk(b,kth1,kth2,mu1,mu2);
    [det_c] = double_fk(c,kth1,kth2,mu1,mu2);
    if ikh > 1 
       if det_b < det_a && det_b < det_c
           [xmin,error] = brent_fk_double(a,b,c,tol,kth1,kth2,mu1,mu2);
           if error < -5
              sum= sum +1;
              root(1,sum)=xmin;
              root_d(1,sum)=error;
           end
       end
    end   
end 

for n_r = 1:size(root,2)
    if root(1,n_r) ~= 0
       kth1 = 2*pi*f/CT1*he;
       kth2 = 2*pi*f/CT2*he;
       [A,D] = double_fk_amplitude(root(1,n_r),kth1,kth2,mu1,mu2 );
       Amp(1:4,n_r) = A;
    else
       break
    end
end

zamp_mode = Amp;
for i = 1:size(root,2)
    kth1 = 2*pi*f/CT1*he;
    kth2 = 2*pi*f/CT2*he;
    kih = root(1,i);
    n_layer = 2;
    miu = [mu1 mu2];
    [ Pm ] = power_eval_t_conju( i,i,n_layer,miu,kth1,kth2,kih,kih,zamp_mode,0,he,1 );
    pm(i,1)=Pm;
    Amp(:,i) = Amp(:,i)/sqrt(abs(Pm));
end

end

