cR = full(cRA);
cT = full(cTS);

figure()
plot(ff,cR(1,:),'-o')
hold on
plot(ff,cR(2,:),'-o')
legend('A','S')

hold off
figure()
plot(ff,cT(1,:),'-or')
hold on
plot(ff,cT(2,:),'-ob')
hold off
legend('A','S')