function [Ke,Me]=Ksolve(E,nu,x1,y1,x2,y2,x3,y3,x4,y4)
%=========单元刚度矩阵===============
% E 弹性模量
% v 泊松比
% h 厚度
% x1,y1,x2,y2,x3,y3,x4,y4 为4个角结点的坐标
%矩阵尺寸为16 x 16
% h=1;
Ke=zeros(16,16);
Me=zeros(16,16);
%% 弹性矩阵
% 平面应变
D = E/(1+nu)/(1-2*nu) * [
        [1-nu nu 0]
        [nu 1-nu 0]
        [0 0 (1- 2*nu)/2]
        ];
%{平面应力
% D = E/(1-nu^2) * [
%         [1 nu 0]
%         [nu 1 0]
%         [0 0 (1- nu)/2]
%         ];
%}

%高斯积分 采用 3 x 3 个积分点 书74页
W1=5/9;W2=8/9;W3=5/9; %加权系数
W=[W1 W2 W3];
r=sqrt(3/5);
x=[-r 0 r];%积分点
for i=1:3
	for j=1:3
		B = eleB(x1,y1,x2,y2,x3,y3,x4,y4,x(i),x(j));
		J = Jacobi(x1,y1,x2,y2,x3,y3,x4,y4,x(i),x(j));
		Ke = Ke+W(i)*W(j)*B'*D*B*det(J);
		Ne = shape(x(i),x(j));
		N = [Ne(1) 0 Ne(2) 0 Ne(3) 0 Ne(4) 0 Ne(5) 0 Ne(6) 0 Ne(7) 0 Ne(8) 0;
			0 Ne(1) 0 Ne(2) 0 Ne(3) 0 Ne(4) 0 Ne(5) 0 Ne(6) 0 Ne(7) 0 Ne(8)];
		Me = Me+W(i)*W(j)*(N')*N*det(J);
	end
end