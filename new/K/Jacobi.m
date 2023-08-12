function J=Jacobi(x1,y1,x2,y2,x3,y3,x4,y4,s,t)

%-------Jacobi-----------
%单元坐标
%5,6,7,8点的坐标
x5=(x1+x2)/2;y5=(y1+y2)/2;
x6=(x2+x3)/2;y6=(y2+y3)/2;
x7=(x3+x4)/2;y7=(y3+y4)/2;
x8=(x1+x4)/2;y8=(y1+y4)/2;
x=[x1 x2 x3 x4 x5 x6 x7 x8];
y=[y1 y2 y3 y4 y5 y6 y7 y8];
%%调用形函数对局部坐标的导数
[N_s,N_t]=DHS(s,t);
%求Jacobi矩阵的行列式的值
x_s=0;y_s=0;
x_t=0;y_t=0;
for i=1:8
x_s=x_s+N_s(i)*x(i);y_s=y_s+N_s(i)*y(i);
x_t=x_t+N_t(i)*x(i);y_t=y_t+N_t(i)*y(i);
end
J=[x_s y_s;x_t y_t];