function plotWave(x,Coordinate,u)
X = Coordinate(:,1);
bnd_node = find(X==x);
bnd = sort([2*bnd_node-1;2*bnd_node]);
uB = u(bnd);
U = uB(1:2:end);
V = uB(2:2:end);
py = Coordinate(bnd_node,2);
[py1,I] = sort(py);
U = U(I);
V = V(I);

figure(Position=[263,382,1106,420])
subplot(1,2,1)
plot(real(U),py1,'b-',LineWidth=1.5);
hold on
plot(imag(U),py1,'b-.')
hold off
subtitle('u1')

subplot(1,2,2)
plot(real(V),py1,'b-',LineWidth=1.5);
hold on
plot(imag(V),py1,'b-.')
hold off
subtitle('u2')
end
