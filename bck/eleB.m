function B=eleB(x1,y1,x2,y2,x3,y3,x4,y4,s,t)
%调用导函数
[N_s,N_t]=DHS(s,t);
%求Jacobi矩阵
J=Jacobi(x1,y1,x2,y2,x3,y3,x4,y4,s,t);
%求应变矩阵B
B=zeros(3,16);
for i=1:8
B1i=J(2,2)*N_s(i)-J(1,2)*N_t(i);
B2i=-J(2,1)*N_s(i)+J(1,1)*N_t(i);
B(1:3,2*i-1:2*i)=[B1i 0;0 B2i;B2i B1i];
end
B=B/det(J);