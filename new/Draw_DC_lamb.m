clear;close all
addpath(".\wave_function")
addpath(".\get_wavenumber")
%% 材料参数
%%% 铝材
lambda = 51e9;
mu = 26e9;
density = 2700;

CL = sqrt((lambda+2*mu)/density);
CT = sqrt(mu/density); 

%% 几何参数
h = 0.5e-3; % 半板厚

%% 扫频范围
f = [.1e6, 5e6];
num_dw = 500; % 扫频点数
w = 2*pi*f;
dw = diff(w)/num_dw; % 扫频步长
w_sca = w(1):dw:w(2) ;

%% 求解波数
Ta = tic; % 计时开始

wd_sca = w_sca*h/CT; % 无量纲化
F1 = @lamb;
[Kd,Wd,nmode] = get_wavenumber(wd_sca,lambda,mu,density,h,F1);

% 量纲还原
K = Kd/h;
W = Wd*CT/h;
F = W/2/pi;

% Cp = W./K;
% Cg = Cp + K.*[nan(1,size(K,2));diff(Cp)./diff(K)];

T = toc(Ta); % 计时结束

%% 绘图
%%% 频散曲线
fig(1) = figure(1);
plot(F(:),K(:),'.');
xlabel("频率f [Hz]")
ylabel("波数k [m^{-1}]")

%%% 相速度
% fig(2) = figure(2);
% plot(F(:),Cp(:),'.')
% xlabel("频率f [Hz]")
% ylabel("相速度C_p [m/s]")
% ylim([0 12e3])

% %%% 群速度
% fig(3) = figure(3);
% plot(F(:),Cg(:),'.')

% xlabel("频率f [Hz]")
% ylabel("群速度C_g [m/s]")
% ylim([0 7e3])

% 保存绘图
% mkdir output % 创建文件夹output，如果文件夹已存在，会有警告，但不影响运行
% date = datetime("now","Format","uuuu-MM-dd HH.mm.ss");
% name_fig = ['./output/Lamb_DC_lamb [', char(date), '].fig'] % 绘图文件名
% savefig(fig,name_fig,'compact')
