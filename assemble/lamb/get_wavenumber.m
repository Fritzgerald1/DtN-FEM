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
Ay = 1e-3;
By = 10;
dy = 1e-3;

%% 搜索波数解
for ii = 1:numel(w_sca) % 扫描频率的区间
	aw = w_sca(ii);
	for  ay = Ay:5*dy:By % 搜索波数解实部的范围
		s1=inf;
		sw1=0;
		sy1=0;
		by = ay+6*dy; % 搜索范围的右边界
		%% 在小网格里，找到极小值解(syr1,siy1,sw1)
		for y1 = ay:dy:by
			h1 = Fun(y1,aw); % 频散方程行列式
			hh1 = abs(h1);
			if hh1 < s1
				s1 = hh1;
				sw1 = aw;
				sy1 = y1;
			end
		end

		if (abs(sy1-ay)>0.1*dy) && (abs(sy1-by)>0.1*dy) % 不在外边界上
			%% 当极小值不在边界时，判断该点是否是零点
			[true_sy,zero_flag]=Determine_zero_point(sy1,sw1,dy,Fun); % 判断零点，ssy1为每次迭代的解
			if zero_flag == 1 % 该点是零点
				Omega=[Omega sw1];
				RealK=[RealK true_sy];
				t=t+1;
			end
		end
	end
end
end