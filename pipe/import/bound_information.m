function [left, right, top ,bottom, freeBnd, left_e, right_e, freeBnd_e] = bound_information(Coordinate, Ielement)
x = Coordinate(:,1);
y = Coordinate(:,2);
L=15*10^(-3);
W=1*10^(-3);

% L_crack = W;

%% �߽�ڵ�
left = find(x==0); % ��߽�ڵ�
right = find(y==min(y)); % �ұ߽�ڵ�
%%	����Ӧ���ڵ�����/��/ȱ�ݱ߽�ڵ����
top = find(y==max(y) .* ~(x==0) .* ~(x==max(x))); % �ϱ߽�ڵ�
% top = find(y==max(y)); % �ϱ߽�ڵ�
% crack = find((y==W/2) .* (x >= L/2-L_crack/2).*(x <= L/2+L_crack/2)); % ȱ�ݱ߽�ڵ�
crack = [];
bottom = find(y==min(y)); % �±߽�ڵ�
freeBnd = [top;crack;bottom]; % ����Ӧ���ڵ�
% freeBnd = [top; bottom];

%% �߽絥Ԫ
[A,~] = arrayfun(@(x) find(ismember(Ielement,x)),left(:),'UniformOutput',false);
left_e = unique(cell2mat(A)); % ��߽絥Ԫ
[A,~] = arrayfun(@(x) find(ismember(Ielement,x)),right(:),'UniformOutput',false);
right_e = unique(cell2mat(A)); % �ұ߽絥Ԫ
%	����Ӧ����Ԫ����/��/ȱ�ݱ߽絥Ԫ���
[A,~] = arrayfun(@(x) find(ismember(Ielement,x)),top(:),'UniformOutput',false);
top_e = unique(cell2mat(A)); % �ϱ߽絥Ԫ
[A,~] = arrayfun(@(x) find(ismember(Ielement(:,5:8),x)),crack(:),'UniformOutput',false);
crack_e = unique(cell2mat(A)); % ȱ�ݱ߽絥Ԫ��	PS��ֻ���ڽڵ�(5~8��)��ȱ����ʱ���õ�Ԫ����ȱ�ݵ�Ԫ
[A,~] = arrayfun(@(x) find(ismember(Ielement,x)),bottom(:),'UniformOutput',false);
bottom_e = unique(cell2mat(A)); % �±߽絥Ԫ
freeBnd_e = [top_e; crack_e; bottom_e]; % ����Ӧ����Ԫ
end