function Amp = get_amplitude(root,wd,lambda,mu,density,h,Fun)
Amp = zeros(4,length(root));

for ii = 1:length(root)
    if root(1,ii) ~= 0
	   [~,~,A] = Fun(root(ii),wd,lambda,mu,density,h);
       Amp(:,ii) = A;
    else
       break
    end
end

end