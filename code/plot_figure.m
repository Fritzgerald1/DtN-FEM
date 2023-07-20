clc;
clear;
f = 1.6:0.1:5.0;
a = load('s_s_Re_m2.mat');
a_name = fieldnames(a);
s_s_Re_m2 = getfield(a,a_name{1});
%value1 = spcrv([[f(1) f f(end)];[s_s_Re_m2(2,1) s_s_Re_m2(2,:) s_s_Re_m2(2,end)]],4);

a = load('s_s_Tr_m2.mat');
a_name = fieldnames(a);
s_s_Tr_m2 = getfield(a,a_name{1});

a = load('s_s_eng_R_m2.mat');
a_name = fieldnames(a);
s_s_eng_R_m2 = getfield(a,a_name{1});

a = load('s_s_eng_R_m1.mat');
a_name = fieldnames(a);
s_s_eng_R_m1 = getfield(a,a_name{1});

a = load('s_s_eng_T_m2.mat');
a_name = fieldnames(a);
s_s_eng_T_m2 = getfield(a,a_name{1});

a = load('s_s_eng_T_m1.mat');
a_name = fieldnames(a);
s_s_eng_T_m1 = getfield(a,a_name{1});

a = load('a_s_Re_m2.mat');
a_name = fieldnames(a);
a_s_Re_m2 = getfield(a,a_name{1});
%value2 = spcrv([[f(1) f f(end)];[a_s_Re_m2(2,1) a_s_Re_m2(2,:) a_s_Re_m2(2,end)]],4);

a = load('a_s_Tr_m2.mat');
a_name = fieldnames(a);
a_s_Tr_m2 = getfield(a,a_name{1});

a = load('ti_s_Re_m2.mat');
a_name = fieldnames(a);
ti_s_Re_m2 = getfield(a,a_name{1});

a = load('ti_s_Tr_m2.mat');
a_name = fieldnames(a);
ti_s_Tr_m2 = getfield(a,a_name{1});


figure
plot(f,s_s_Re_m2(1,:));
hold on;
plot(f,s_s_Re_m2(2,:),'b--o');
hold on;
%plot(value1(1,:),value1(2,:));
hold on;
plot(f,a_s_Re_m2(1,:));
hold on;
plot(f,a_s_Re_m2(2,:),'g--x');
%plot(value2(1,:),value2(2,:));
hold on;
plot(f,ti_s_Re_m2(1,:));
hold on;
plot(f,ti_s_Re_m2(2,:),'r--*');

figure
plot(f,s_s_Tr_m2(1,:));
hold on;
plot(f,s_s_Tr_m2(2,:));
hold on;
plot(f,a_s_Tr_m2(1,:));
hold on;
plot(f,a_s_Tr_m2(2,:));
hold on;
plot(f,ti_s_Tr_m2(1,:));
hold on;
plot(f,ti_s_Tr_m2(2,:));

figure
plot(f,s_s_eng_R_m1(1,:));
hold on;
plot(f,s_s_eng_T_m1(1,:));
hold on;
plot(f,s_s_eng_R_m2(1,:));
hold on;
plot(f,s_s_eng_T_m2(1,:));