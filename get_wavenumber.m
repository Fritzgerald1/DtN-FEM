function [RealK,Omega,t] = get_wavenumber(w_sca,Fun)
%%
% 输入量：
%	w_sca为无量纲数
%		w_sca = w*h/CT
% 输出量：
%	RealK也是无量纲数
%		RealK = k*h
%	t是模态数量
%	Omega是t个w_sca

%% 变量初始化
RealK = [];
Omega = [];
t=0;

%% 搜索波数解的范围
Ary = -1e-2;
Bry = 7;
dry = 1e-2;


%% 搜索波数解
parfor ii = 1:numel(w_sca) % 扫描频率的区间
	aw = w_sca(ii);
	for  ary = Ary:10*dry:Bry % 搜索波数解实部的范围
		s1=inf;
		sw=0;
		sry=0;
		bry = ary+11*dry; % 搜索范围的实部右边界
		%% 在小网格里，找到极小值解(syr1,siy1,sw1)
		for ry1 = ary:dry:bry
			y1 = ry1;
			h1 = Fun(y1,aw); % 频散方程行列式
			hh1 = abs(h1);
			if hh1 < s1
				s1 = hh1;
				sw = aw;
				sry = ry1;
			end
		end

		if (abs(sry-ary)>0.1*dry) && (abs(sry-bry)>0.1*dry) % 不在外边界上
			%% 当极小值不在边界时，判断该点是否是零点
			[zsy,zero_flag]=Determine_zero_point(sry,sw,dry,Fun); % 判断零点，ssy1为每次迭代的解
			if zero_flag == 1 % 该点是零点
				Omega=[Omega sw];
				RealK=[RealK zsy];
				t=t+1;
			end
		end
	end
end