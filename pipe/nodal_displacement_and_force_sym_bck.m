function [U,T] = nodal_displacement_and_force(wd,Kd,n,Amp,bnd,Coordinate,lambda,mu,density,h, alpha)
%% ���ұ߽紦�Ľڵ�ģ̬�ڵ�λ�Ʀ�����T (��n��ģ̬)
CL = sqrt((lambda+2*mu)/density);
CT = sqrt(mu/density); 

kd = Kd(n); % ��k�ײ���
k = kd/h;
w = wd*CT/h;

p = sqrt((w/CL)^2-k.^2);
q = sqrt((w/CT)^2-k.^2);

X = Coordinate(bnd,1);
Y = Coordinate(bnd,2);


E =@(x,y) [exp(1i*p*y) exp(-1i*p*y) exp(1i*q*y) exp(-1i*q*y)];
A = Amp(:,n); % ��n��ģ̬�ķ�ֵϵ��

u = @(x,y) E(x,y) * diag([k, k, -q, q].*1i) * A .* exp(1i*k*x);
v = @(x,y) E(x,y) * diag([p, -p k k].*1i) * A .* exp(1i*k*x);
u = @(x,y) imag(u(x,y)); v =@(x,y) real(v(x,y));

Ux =@(x,y) cos(alpha)*u(x,y) + sin(alpha)*v(x,y);
Uy =@(x,y) sin(alpha)*u(x,y) + cos(alpha)*v(x,y);

t11 =@(x,y) E(x,y) * diag([mu*(2*p^2-q^2-k^2)*[1 1], 2*mu*q*k*[1 -1]]) * A .* exp(1i*k*x);
t12 =@(x,y) E(x,y) * diag([2*mu*k*p*[-1 1], mu*(q^2-k^2)*[1 1]]) * A .* exp(1i*k*x);
t22 =@(x,y) E(x,y) * diag([-mu*(q^2-k^2)*[1 1], 2*mu*k*q*[-1 1]]) * A .* exp(1i*k*x);
t11 =@(x,y) real(t11(x,y)); t12 =@(x,y) imag(t12(x,y)); t22 =@(x,y) real(t22(x,y)); 
tx = @(x,y) t11(x,y)+t12(x,y);
ty = @(x,y) t12(x,y)+t22(x,y);
Tx = @(x,y) cos(alpha)*tx(x,y) + sin(alpha)*ty(x,y);
Ty = @(x,y) sin(alpha)*tx(x,y) + cos(alpha)*ty(x,y);

% �ϲ�x��y��λ��/��
Ux1 = Ux(X,Y);
Uy1 = Uy(X,Y);
tmpU = [Ux1,Uy1]';
U = tmpU(:);

Tx1 = Tx(X,Y);
Ty1 = Ty(X,Y);
tmpT = [Tx1,Ty1]';
T = tmpT(:);



% %% ��ͼ����
%λ��
figure;
[y1,I] = sort(Y);
plot(Ux1(I),y1)
grid on
hold on
plot(Uy1(I),y1,'--')
legend('u','v')
%��
figure;
[y1,I] = sort(Y);
plot(Tx1(I),y1)
grid on
hold on
plot(Ty1(I),y1,'--')
legend('tx','ty')

end