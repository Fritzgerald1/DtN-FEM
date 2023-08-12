%%

dof = 2*Nnode;
u = ones(dof,1);

%% 定义边界
% left = [1 4 8]; % 边界节点编号
% right = [2 3 6];
% bottom = [1 2 5];
% top = [3 4 7];

left = bnd_L;
right = bnd_R;
top = bnd_T;
bottom = bnd_B;

lx = 2*left-1; ly = 2*left; % 边界x、y方向的行号
rx = 2*right-1; ry = 2*right;
bx = 2*bottom-1; by = 2*bottom;
tx = 2*top-1; ty = 2*top;

%% 初始位移
u(lx)=0; % 左端u=0
% u(10) = 0; % 节点5处v=0
u(by)=0; % 底端v=0

%% 边界条件

P = [1e10 0]'; % 右端外力

% 等效节点力
F = zeros(dof,1);
w=[1 1]; r=.577; % 加权系数和位置
x=[-r r]; % 积分点
for i=1:2
	Ne=shape(1,x(i));
	N = [Ne(1) 0 Ne(2) 0 Ne(3) 0 Ne(4) 0 Ne(5) 0 Ne(6) 0 Ne(7) 0 Ne(8) 0;
		0 Ne(1) 0 Ne(2) 0 Ne(3) 0 Ne(4) 0 Ne(5) 0 Ne(6) 0 Ne(7) 0 Ne(8)];
	F = F + w(i)*N'*P*(1/2);
end


for i = 1:8
	right_e([2*i-1 2*i]) = [2*Ielement(bnd_R_e,i)-1,2*Ielement(bnd_R_e,i)];
end
T(right_e,1) = F(:);
%

%% 求解未知自由度
index=find(u); %未知自由度的索引
% index=[];
% for i=1:dof
% 	if u(i)~=0
% 		index=[index,i];
% 	end
% end
u(index) = Kg(index,index)\T(index);

%% 变形图
figure()
ux = u(1:2:end);
uy = u(2:2:end);

% px = [0 1 1 0 .5 1 .5 0]';
% py = [0 0 1 1 0 .5 1 .5]';
px = [0 0 1 1 0 .5 .5 1]';
py = [-.5 .5 -.5 .5 0 -.5 .5 0]';

dx = px+ux;
dy = py+uy;
% order = [1 5 2 6 3 7 4 8 1];
order = [1 6 3 8 4 7 2 5 1];
plot(px(order),py(order),'k')
% plot(px(:),py(:),'k')
hold on
plot(dx(order),dy(order),'r')
% plot(dx(:),dy(:),'r')
for i=1:9
	text(dx(order(i)),dy(order(i)),num2str(order(i)))
end
hold off
axis padded


%% 查看右端x方向位移u
figure()
ru = u(rx);
rv = u(ry);
plot(ru,1:3)
hold on
plot(rv,1:3)
hold off
legend('u','v')
yticks([1 2 3]); grid on