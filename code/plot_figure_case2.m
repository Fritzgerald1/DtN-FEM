clc;
clear;
f = 1.6:0.1:5.0;
a = load('s_s_Re_m2.mat');
a_name = fieldnames(a);
s_s_Re_m2 = getfield(a,a_name{1});

a = load('a_s_Re_m2.mat');
a_name = fieldnames(a);
a_s_Re_m2 = getfield(a,a_name{1});

a = load('ti_s_Re_m2.mat');
a_name = fieldnames(a);
ti_s_Re_m2 = getfield(a,a_name{1});

a = load('s_s_Tr_m2.mat');
a_name = fieldnames(a);
s_s_Tr_m2 = getfield(a,a_name{1});

a = load('a_s_Tr_m2.mat');
a_name = fieldnames(a);
a_s_Tr_m2 = getfield(a,a_name{1});

a = load('ti_s_Tr_m2.mat');
a_name = fieldnames(a);
ti_s_Tr_m2 = getfield(a,a_name{1});

a = load('s_s_eng_R_m2.mat');
a_name = fieldnames(a);
s_s_eng_R_m2 = getfield(a,a_name{1});

a = load('a_s_eng_R_m2.mat');
a_name = fieldnames(a);
a_s_eng_R_m2 = getfield(a,a_name{1});

a = load('ti_s_eng_R_m2.mat');
a_name = fieldnames(a);
ti_s_eng_R_m2 = getfield(a,a_name{1});

a = load('s_s_eng_T_m2.mat');
a_name = fieldnames(a);
s_s_eng_T_m2 = getfield(a,a_name{1});

a = load('a_s_eng_T_m2.mat');
a_name = fieldnames(a);
a_s_eng_T_m2 = getfield(a,a_name{1});

a = load('ti_s_eng_T_m2.mat');
a_name = fieldnames(a);
ti_s_eng_T_m2 = getfield(a,a_name{1});

figure(1)
plot(f,s_s_Re_m2(1,:),f,a_s_Re_m2(1,:),f,ti_s_Re_m2(1,:));
hold on;
plot(f,s_s_Re_m2(2,:),f,a_s_Re_m2(2,:),f,ti_s_Re_m2(2,:));
legend('s-s-inv2-m1','a-s-inv2-m1','ti-s-inv2-m1','s-s-inv2-m2','a-s-inv2-m2','ti-s-inv2-m2');

figure(2)
plot(f,s_s_Tr_m2(1,:),f,a_s_Tr_m2(1,:),f,ti_s_Tr_m2(1,:));
hold on;
plot(f,s_s_Tr_m2(2,:),f,a_s_Tr_m2(2,:),f,ti_s_Tr_m2(2,:));
legend('s-s-inv2-m1','a-s-inv2-m1','ti-s-inv2-m1','s-s-inv2-m2','a-s-inv2-m2','ti-s-inv2-m2');

figure(3)
plot(f,s_s_eng_R_m2,f,a_s_eng_R_m2,f,ti_s_eng_R_m2);
hold on;
plot(f,s_s_eng_T_m2,f,a_s_eng_T_m2,f,ti_s_eng_T_m2);
legend('s-s-inv2-engR','a-s-inv2-engR','ti-s-inv2-engR','s-s-inv2-engT','a-s-inv2-engT','ti-s-inv2-engT');