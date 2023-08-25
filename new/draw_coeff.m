close all
addpath("./output")
cr = abs(cR);
ct = abs(cT);
h = 1e-3;
ff1 = ff/CT*h*2*pi/pi*2;

fig(1) = figure(Name='1',Position=[127.4,221,1279.2,420]);
subplot(1,2,1)
plot(ff1,cr(1,:),'-o','Marker','.')
hold on
plot(ff1,cr(2,:),'-o','Marker','.')
hold off
legend('A0','S0')
title('Reflection')
axis padded

% h(2) = figure(2);
subplot(1,2,2)
plot(ff1,ct(1,:),'-o','Marker','.')
hold on
plot(ff1,ct(2,:),'-o','Marker','.')
hold off
legend('A0','S0')
axis padded
title('Transmission')

char1 = 'A0'; % 额外的命名字段
date = datetime("now","Format","uuuu-MM-dd HH.mm.ss");
name_fig = ['./output/',char1 , '[', char(date), '].fig'] % 绘图文件名
savefig(fig,name_fig)