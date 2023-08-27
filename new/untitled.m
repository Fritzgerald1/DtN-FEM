x=0;
X = Coordinate(:,1);
bnd = find(X==x);

b = sort([2*bnd-1;2*bnd]);
y = Coordinate(bnd,2);


u = nan(2*Nnode,1);
u(b) = u;

uB = u(b);
U = uB(1:2:end);
V = uB(2:2:end);

[yi,I] = sort(y);
U = U(I);
V = V(I);

figure
plot(real(U),yi,'b-',LineWidth=1.5);
hold on
plot(real(V),yi,'r-.')
hold off
% 
% figure
% plot(real(V),yi,'b-',LineWidth=1.5);
% hold on
% plot(imag(V),yi,'r-.')
% hold off
