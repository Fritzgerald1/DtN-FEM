function [left, right, top ,bottom, freeBnd, left_e, right_e, freeBnd_e] = bound_information(Coordinate, Ielement)
x = Coordinate(:,1);
y = Coordinate(:,2);
L=15*10^(-3);
W=1*10^(-3);

% L_crack = W;

%% 边界节点
left = find(x==0); % 左边界节点
right = find(y==min(y)); % 右边界节点
%%	自由应力节点由上/下/缺陷边界节点组成
top = find(y==max(y) .* ~(x==0) .* ~(x==max(x))); % 上边界节点
% top = find(y==max(y)); % 上边界节点
% crack = find((y==W/2) .* (x >= L/2-L_crack/2).*(x <= L/2+L_crack/2)); % 缺陷边界节点
crack = [];
bottom = find(y==min(y)); % 下边界节点
freeBnd = [top;crack;bottom]; % 自由应力节点
% freeBnd = [top; bottom];

%% 边界单元
[A,~] = arrayfun(@(x) find(ismember(Ielement,x)),left(:),'UniformOutput',false);
left_e = unique(cell2mat(A)); % 左边界单元
[A,~] = arrayfun(@(x) find(ismember(Ielement,x)),right(:),'UniformOutput',false);
right_e = unique(cell2mat(A)); % 右边界单元
%	自由应力单元由上/下/缺陷边界单元组成
[A,~] = arrayfun(@(x) find(ismember(Ielement,x)),top(:),'UniformOutput',false);
top_e = unique(cell2mat(A)); % 上边界单元
[A,~] = arrayfun(@(x) find(ismember(Ielement(:,5:8),x)),crack(:),'UniformOutput',false);
crack_e = unique(cell2mat(A)); % 缺陷边界单元。	PS：只有内节点(5~8号)在缺陷上时，该单元才算缺陷单元
[A,~] = arrayfun(@(x) find(ismember(Ielement,x)),bottom(:),'UniformOutput',false);
bottom_e = unique(cell2mat(A)); % 下边界单元
freeBnd_e = [top_e; crack_e; bottom_e]; % 自由应力单元
end