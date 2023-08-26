function [Ko,Kl,Kr,l,r,o] = Stiffness_Changed_Position(Kg,bnd_L,bnd_R,Nnode)
%% ��������նȾ���Kg��ʹ��λ�ƺ�Ӧ��˳��Ϊ[u0;uL;uR]��[f0;tL;tR]

%{
% �����Ƿ񽻻�׼ȷ
Kg = cell(2*Nnode);
for i = 1:2*Nnode
	for j =1:2*Nnode
		Kg{i,j}=[i,j];
	end
end
%}

l = sort([2*bnd_L-1;2*bnd_L]); % ��߽�ڵ���к�
r = sort([2*bnd_R-1;2*bnd_R]); % �ұ߽�ڵ���к�
o = 1:2*Nnode; 
o([l;r])=[]; % �м�߽�ڵ���к�

l = l'; r = r'; o = o';
Ko = [Kg(o,o);Kg(l,o);Kg(r,o)];
Kl = [Kg(o,l);Kg(l,l);Kg(r,l)];
Kr = [Kg(o,r);Kg(l,r);Kg(r,r)];
end