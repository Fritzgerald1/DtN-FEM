function [Ko,Kl,Kr,l,r,o] = Stiffness_Changed_Position(Kg,L,R,Nnode)
%% 重新排序刚度矩阵Kg，使得位移和应力顺序为[u0;uL;uR]和[f0;tL;tR]

%{ 测试是否交换准确
% Kg = cell(2*Nnode);
% for i = 1:2*Nnode
% 	for j =1:2*Nnode
% 		Kg{i,j}=[i,j];
% 	end
% end
%}

% L = bnd_401; % 左边界节点号
% R = bnd_501; % 右边界节点号
l = sort([2*L-1;2*L]); % 左边界节点的行号
r = sort([2*R-1;2*R]); % 右边界节点的行号
o = 1:2*Nnode; 
o([l;r])=[]; % 非左右边界节点的行号

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