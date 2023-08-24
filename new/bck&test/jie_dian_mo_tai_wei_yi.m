% ���䲨�ұ߽��ģ̬�ڵ�λ�Ʀ���ģ̬��T(����ڶ���ģ̬,Ƶ��wd,����kd,)

%% ������
n = 2; % ģ̬����
k = K(n);
CL = 6.35e3;
CT = 3.13e3;
mu = 26e9;

p = sqrt((w/CL)^2-k^2);
q = sqrt((w/CT)^2-k^2);

% x = linspace(0,15e-3,1000);
x = Coordinate(bnd_R,1);
y = Coordinate(bnd_R,2);
% y = y(1);

E = [exp(1i*p*y) exp(-1i*p*y) exp(1i*q*y) exp(-1i*q*y)];
A = Amp(:,n); % ��n��ģ̬�ķ�ֵϵ��

u =  E * diag([k, k, -q, q].*1i) * A .* exp(1i*k*x);
v =  E * diag([p, -p k k].*1i) * A .* exp(1i*k*x);
t11 = E * diag([mu*(2*p^2-q^2-k^2)*[1 1], 2*mu*q*k*[1 -1]]) * A .* exp(1i*k*x);
t12 = E * diag([2*mu*k*p*[-1 1], mu*(q^2-k^2)*[1 1]]) * A .* exp(1i*k*x);
t22 = E * diag([-mu*(q^2-k^2)*[1 1], 2*mu*k*q*[-1 1]]) * A .* exp(1i*k*x);

% %% ��ͼ����
% figure
% [y1,I] = sort(y);
% plot(u(I),y1)
% % plot(x,v)
end