function Ke = Elem801(lambda,mu,h,x_e,y_e)

E = mu*(2*mu+3*lambda)/(mu+lambda);
nu = lambda/(2*(mu+lambda));

    syms x y
    % 对于平面八节点单元，共有
    M = [
        [1 x y x^2 x*y y^2 x^2*y x*y^2 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 1 x y x^2 x*y y^2 x^2*y x*y^2]
        ];
    
    
    % 定义八节点单元的x_e和y_e
    % x_e = [-1 1 1 -1 0 1 0 -1];
    % y_e = [-1 -1 1 1 -1 0 1 0];
    
    % 创建初始矩阵A
    A = zeros(16,16);
    for i = 1:8
        A(2*i-1:1:2 *i,:) = subs(M,[x,y],[x_e(i),y_e(i)]);
    end
    
    % 计算形状函数矩阵 
    N = M/A;
    
    % 使用求导方法求解对应的位移函数
    N_x = diff(N,x);
    N_y = diff(N,y);
    B = [N_x(1,:)
         N_y(2,:)
         N_y(1,:) + N_x(2,:)
         ];
    clear N_x N_y
    % 刚度矩阵D
    D = E/(1-nu^2) * [
        [1 nu 0]
        [nu 1 0]
        [0 0 (1- nu)/2]
        ];


    
    ke_pre = transpose(B)*D*B;  % 注意不能使用B'
    Ke = h * int(int(ke_pre,x,[-1, 1]), y, [-1,1]); % 积分得到刚度矩阵
    % 注意不能使用vpa函数, vpa函数会产生截断误差
	Ke = double(Ke);
    clear ke_pre
end