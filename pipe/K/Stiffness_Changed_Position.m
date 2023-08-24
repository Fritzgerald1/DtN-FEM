function [Ko,Kl,Kr,l,r,o] = Stiffness_Changed_Position(Kg,L,R,Nnode)
%% ��������նȾ���Kg��ʹ��λ�ƺ�Ӧ��˳��Ϊ[u0;uL;uR]��[f0;tL;tR]

%{ �����Ƿ񽻻�׼ȷ
% Kg = cell(2*Nnode);
% for i = 1:2*Nnode
% 	for j =1:2*Nnode
% 		Kg{i,j}=[i,j];
% 	end
% end
%}

% L = bnd_401; % ��߽�ڵ��
% R = bnd_501; % �ұ߽�ڵ��
l = sort([2*L-1;2*L]); % ��߽�ڵ���к�
r = sort([2*R-1;2*R]); % �ұ߽�ڵ���к�
o = 1:2*Nnode; 
o([l;r])=[]; % �����ұ߽�ڵ���к�

Kg_oo = Kg(o,o);
Kg_ol = Kg(o,l);
Kg_or = Kg(o,r);
Kg_ll = Kg(l,l);
Kg_rr = Kg(r,r);
Kg_lr = Kg(l,r);

% Kg_new = [Kg_oo,Kg_ol,Kg_or;
% 		  Kg_ol',Kg_ll,Kg_lr;
% 		  Kg_or',Kg_lr',Kg_rr];
Ko = [Kg_oo;Kg_ol';Kg_or'];
Kl = [Kg_ol;Kg_ll;Kg_lr'];
Kr = [Kg_or;Kg_lr;Kg_rr];
end