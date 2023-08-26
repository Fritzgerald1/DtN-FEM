%% 重新排序刚度矩阵Kg，使得位移和应力顺序为[u0;uL;uR]和[f0;tL;tR]

% 测试是否交换准确
Kg = cell(2*Nnode);
for i = 1:2*Nnode
	for j =1:2*Nnode
		Kg{i,j}=[i,j];
	end
end
%}
L = bnd_L;
R = bnd_R;
% L = bnd_401; % 左边界节点号
% R = bnd_501; % 右边界节点号
l = sort([2*L-1;2*L]); % 左边界节点的行号
r = sort([2*R-1;2*R]); % 右边界节点的行号
o = 1:2*Nnode; 
o([l;r])=[]; % 非左右边界节点的行号

l = l';
r = r';
o = o';
Ko = [Kg(o,o);Kg(l,o);Kg(r,o)];
Kl = [Kg(o,l);Kg(l,l);Kg(r,l)];
Kr = [Kg(o,r);Kg(l,r);Kg(r,r)];