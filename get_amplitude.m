function [Amp zamp_mode] = get_amplitude(root,wd,he,Fun)
Amp = zeros(4,length(root));

for ii = 1:length(root)
    if root(1,ii) ~= 0
	   [~,~,A] = Fun(root(ii),wd);
       Amp(:,ii) = A;
    else
       break
    end
end


zamp_mode = Amp;
for i = 1:length(root)
    kih = root(1,i);
    n_layer = 1;
% 	[ Pm ] = power_eval_t_conju( i,i,n_layer,mu,wd,wd,kih,kih,zamp_mode,0,he,1 );
%     [ Pm ] = power_eval_t_conju_new( i,i,n_layer,mu,wd,wd,kih,zamp_mode,0,he,1 );
%     pm(i,1)=Pm;
%     Amp(:,i) = Amp(:,i)/sqrt(abs(Pm));
end
end