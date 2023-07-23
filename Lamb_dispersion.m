function [ zamp_mode,Amp ] = Lamb_dispersion( mu,density,f,he,root,wd)
CT1 = sqrt(mu/density);
Amp = zeros(2,size(root,2));

for ii = 1:size(root,2)
    if root(1,ii) ~= 0
	   [A] = get_amp(root(ii),wd);
%        [A,D] = amplitude(root(1,n_r),kth1,kth1,mu,mu );
       Amp(:,ii) = A;
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
    miu = [mu mu2];
    [ Pm ] = power_eval_t_conju( i,i,n_layer,miu,kth1,kth2,kih,kih,zamp_mode,0,he,1 );
    pm(i,1)=Pm;
    Amp(:,i) = Amp(:,i)/sqrt(abs(Pm));
end

end

