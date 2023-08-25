close all
clear h
addpath("./output")
cr = abs(cR);
ct = abs(cT);
h = 1e-3;
ff1 = ff/CT*h*2*pi/pi*2;

h(1) = figure(1);
plot(ff1,cr(1,:),'-o','Marker','.')
hold on
plot(ff1,cr(2,:),'-o','Marker','.')
legend('A','S')
title('Reflection')
hold off


h(2) = figure(2);
plot(ff1,ct(1,:),'-o','Marker','.')
hold on
plot(ff1,ct(2,:),'-o','Marker','.')
hold off
legend('A','S')
title('Transmission')

char1 = 'crack S0'; % 额外的命名字段
date = datetime("now","Format","uuuu-MM-dd HH.mm.ss");
name_fig = ['./output/',char1 , '[', char(date), '].fig'] % 绘图文件名
savefig(h,name_fig)



% for tt = 1:nff
% f = ff(tt); % 入射频率
% w = 2*pi*f;
% mode_in = 1; % A0模态入射
% wd = w/CT*h; % 无量纲频率
% Fun = @lamb;
% [kd,~,modes] = get_wavenumber(wd,lambda,mu,density,h,Fun); % 得到无量纲波数和模态个数
% kd = kd(end:-1:1); % 波数按低阶到高阶排列
% k = kd/h; % 去归一化
% [Amp] = get_amplitude( kd,wd,lambda,mu,density,h,Fun );
% CR(tt) = sum(abs(Amp(:,1)*cr(1,tt)+Amp(:,2)*cr(2,tt)))./sum(abs(Amp(:,mode_in)));
% end
% 
% plot(ff,CR,'-o','Marker','.')