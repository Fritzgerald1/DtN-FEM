clc;
clear;
f = 1.6:0.1:5.0;
a = load('a_s_all_BEM_600.mat');
a_name = fieldnames(a);
BEM = getfield(a,a_name{1});

a = load('a_s_Re_m1.mat');
a_name = fieldnames(a);
a_s_Re_m1 = getfield(a,a_name{1});

a = load('a_s_Tr_m1.mat');
a_name = fieldnames(a);
a_s_Tr_m1 = getfield(a,a_name{1});

a = load('a_s_eng_R_m1.mat');
a_name = fieldnames(a);
a_s_eng_R_m1 = getfield(a,a_name{1});

a = load('a_s_eng_T_m1.mat');
a_name = fieldnames(a);
a_s_eng_T_m1 = getfield(a,a_name{1});


a = load('a_s_Re_m2.mat');
a_name = fieldnames(a);
a_s_Re_m2 = getfield(a,a_name{1});

a = load('a_s_Tr_m2.mat');
a_name = fieldnames(a);
a_s_Tr_m2 = getfield(a,a_name{1});

a = load('a_s_eng_R_m2.mat');
a_name = fieldnames(a);
a_s_eng_R_m2 = getfield(a,a_name{1});

a = load('a_s_eng_T_m2.mat');
a_name = fieldnames(a);
a_s_eng_T_m2 = getfield(a,a_name{1});

figure(1)
for i = 1:size(BEM,1)-2
    if mod(i,2) ~= 0
        plot(f,BEM(i,:))
        hold on;
    end
end

for i = 1:2
    plot(f,a_s_Re_m1(i,:));
    hold on;
    plot(f,a_s_Re_m2(i,:));
    hold on;
end
legend('BEM-inv2-m1','BEM-inv2-m2','FEM-inv1-m1','FEM-inv2-m1','FEM-inv1-m2','FEM-inv2-m2');
%set(gca,'looseInset',[0 0 0 0]);

figure(2)
for i = 1:size(BEM,1)-2
    if mod(i,2) == 0
        plot(f,BEM(i,:))
        hold on;
    end
end

for i = 1:2
    plot(f,a_s_Tr_m1(i,:));
    hold on;
    plot(f,a_s_Tr_m2(i,:));
    hold on;
end
legend('BEM-inv2-m1','BEM-inv2-m2','FEM-inv1-m1','FEM-inv2-m1','FEM-inv1-m2','FEM-inv2-m2');

figure(3)
for i = size(BEM,1)-1:size(BEM,1)
   plot(f,BEM(i,:))
   hold on;
end
plot(f,a_s_eng_R_m1,f,a_s_eng_T_m1);
hold on;
plot(f,a_s_eng_R_m2,f,a_s_eng_T_m2);
legend('BEM-inv2-engR','BEM-inv2-engT','FEM-inv1-engR','FEM-inv1-engT','FEM-inv2-engR','FEM-inv2-engT');
