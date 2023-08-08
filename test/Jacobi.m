function J=Jacobi(ie,s,t,Elements,Nodes)
ENodes = Elements(ie,:); 
%获取单元结点
xe = Nodes(ENodes(:),:);
%获取节点坐标
x1=xe(1,1);
y1=xe(1,2);
x2=xe(2,1);
y2=xe(2,2);
x3=xe(3,1);
y3=xe(3,2);
x4=xe(4,1);
y4=xe(4,2);
J=1/4*[-(1+t), -(1-t), 1-t, 1+t;
	1-s, -(1-s), -(1+s), 1+s] * [x1 y1;x2 y2;x3 y3;x4 y4];
end