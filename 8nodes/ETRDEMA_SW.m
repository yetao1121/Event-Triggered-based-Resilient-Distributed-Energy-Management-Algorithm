function [SW, k_axis] = ETRDEMA_SW()

clear;clc;
n = 8;% 节点数量

% 发电端参数
a = [0.0024, 0.0056, 0.0072];
b = [5.56, 4.32, 6.60];
c = [30, 25, 25];

% 传输损耗参数：
B = [0.000034 0.000013 0.000009;
     0.000013 0.000014 0.000010;
     0.000009 0.000010 0.000031];
B0 = [-0.00000075 -0.00000006 0.00000070];
B00 = 0.1; 

% 负载端参数
w = [0, 0, 0, 18.43, 13.17, 15.46, 10.03, 8.45];
e = [0, 0, 0, 0.0545, 0.0877, 0.0547, 0.1041, 0.0870];

% 局部约束上下限
p_min = [60, 25, 28, 50, 100, 40, 30, 80];
p_max = [339.69, 479.10, 290.4, 100.34, 159.13, 80.56, 123.98, 109.55];

% 初始化

%初始化每个节点的lambda值
for i = 1 : 1 : n
    if i <= 3
        lambda_state(i) = (2 * a(i) * p_max(i) + b(i));
        lambda_trig(i) = lambda_state(i);
    else
        lambda_state(i) = -2 * e(i) * p_min(i) + w(i);
        lambda_trig(i) = lambda_state(i);
    end
end

% 初始化每个节点的功率
P_state = zeros(1,n);

% 初始化每个节点的 ksai 值
ksai_state = zeros(1,n);
ksai_trig = zeros(1,n);

% 记录每个节点的触发次数
count = zeros(1,n);

% 设定迭代次数
k_max = 800;

% 没有事件触发每个节点都是 k_max 次触发，做柱状图对比
count_notrig = k_max * ones(1, n);

% 循环控制参数设置
eta = 0.000708;

% 定义数组存放每一次迭代的数据以用于作图
for i = 1 : 1 : n
    P_axis{i} = 0;
    Lambda{i} = 0;
    Ksai{i} = 0;
    k_trig{i} = [];
end
k_axis = 0;
ploss = 0;
delta_P = 0;
SW = 0;

for k = 1 : 1 : k_max
    % 迭代步长
    ak = 1 / (k + 1);
    
    % 事件触发阈值
    Ek = 1 / (0.5 * k + 1);
    
    % 每个节点更新 lambda 值
    for i = 1 : 1 : n
        %定义包含自身以及邻居的状态集（其中邻居为上一次触发时刻状态，自身为实时状态）
        Nl{1} = [lambda_state(1),lambda_trig(4),lambda_trig(5)];
        Nl{2} = [lambda_state(2),lambda_trig(4)];
        Nl{3} = [lambda_state(3),lambda_trig(5)];
        Nl{4} = [lambda_trig(1),lambda_trig(2),lambda_state(4),lambda_trig(6),lambda_trig(7)];
        Nl{5} = [lambda_trig(1),lambda_trig(3),lambda_state(5),lambda_trig(8)];
        Nl{6} = [lambda_trig(4),lambda_state(6)];
        Nl{7} = [lambda_trig(4),lambda_state(7)];
        Nl{8} = [lambda_trig(5),lambda_state(8)];

        %定义包含自身以及可信邻居的状态集（其中邻居为上一次触发时刻状态，自身为实时状态）
        Tl{1} = [lambda_state(1),lambda_trig(4),lambda_trig(5)];
        Tl{2} = [lambda_trig(4),lambda_state(2)];
        Tl{3} = [lambda_trig(5),lambda_state(3)];
        Tl{4} = [lambda_trig(1),lambda_state(4)];
        Tl{5} = [lambda_trig(1),lambda_state(5)];
        Tl{6} = [lambda_trig(4),lambda_state(6)];
        Tl{7} = [lambda_trig(4),lambda_state(7)];
        Tl{8} = [lambda_trig(5),lambda_state(8)];
        
        Tl{i} = sort(Tl{i});
        R = [];
        
        % 有效信息筛选
        for j = 1 : 1 : length(Nl{i})
            if Nl{i}(j) <= Tl{i}(end) && Nl{i}(j) >= Tl{i}(1)
                R(end + 1) = Nl{i}(j);
            end
        end
        
        % lambda_i 开始更新
        sigma = 0;
        for j = 1 : 1 : length(R)
            sigma = sigma + R(j);
        end
        sigma = sigma / length(R);
        lambda_state(i) = sigma + eta * ksai_state(i);
        
        % 对 lambda 进行攻击模拟
        lambda_state(6) = 4;
        lambda_state(7) = 4 * sin(0.004 * pi * k) + 4;
        lambda_state(8) = (k/300)^2;
        
        Lambda{i}(end + 1) = lambda_state(i);
    end
    
    % 本地更新功率 P    
    for i = 1 : 1 : n
        if i <= 3
            P_state(i) = motorfunc(i, lambda_state(i), p_min(i), p_max(i));
        else
            P_state(i) = loadfunc(i, k, lambda_state(i), p_min(i), p_max(i));
        end
        P_axis{i}(end + 1) = P_state(i);
    end
    
    % 计算目标函数值，即当前社会福利值
    sw = 0;
    for i = 1 : 1 : 3
        sw = sw + a(i) * P_state(i)^2 + b(i) * P_state(i) + c(i);
    end
    for i = 4 : 1 : 8
        if P_state(i) <= w(i) / (2 * e(i))
            sw = sw - w(i) * P_state(i) + e(i) * P_state(i);
        else
            sw = sw - w(i)^2 / (4 * e(i));
        end
    end
    SW(end + 1) = sw;
    
    % 计算 ploss
    pls = 0;
    for i = 1 : 1 : 3
        for j = 1 : 1 : 3
            pls = pls + B(i,j) * P_state(i) * P_state(j);
        end
    end
    for i = 1:1:3
        pls = pls + B0(i) * P_state(i);
    end
    pls = pls + B00;
    ploss(end + 1) = pls;
    
    % 每个节点更新 ksai 值
    for i = 1 : 1 : n
        %定义包含自身以及邻居的状态集（其中邻居为上一次触发时刻状态，自身为实时状态）
        Nk{1} = [ksai_state(1),ksai_trig(4),ksai_trig(5)];
        Nk{2} = [ksai_state(2),ksai_trig(4)];
        Nk{3} = [ksai_state(3),ksai_trig(5)];
        Nk{4} = [ksai_trig(1),ksai_trig(2),ksai_state(4),ksai_trig(6),ksai_trig(7)];
        Nk{5} = [ksai_trig(1),ksai_trig(3),ksai_state(5),ksai_trig(8)];
        Nk{6} = [ksai_trig(4),ksai_state(6)];
        Nk{7} = [ksai_trig(4),ksai_state(7)];
        Nk{8} = [ksai_trig(5),ksai_state(8)];

        %定义包含自身以及可信邻居的状态集（其中邻居为上一次触发时刻状态，自身为实时状态）
        Tk{1} = [ksai_state(1),ksai_trig(4),ksai_trig(5)];
        Tk{2} = [ksai_trig(4),ksai_state(2)];
        Tk{3} = [ksai_trig(5),ksai_state(3)];
        Tk{4} = [ksai_trig(1),ksai_state(4)];
        Tk{5} = [ksai_trig(1),ksai_state(5)];
        Tk{6} = [ksai_trig(4),ksai_state(6)];
        Tk{7} = [ksai_trig(4),ksai_state(7)];
        Tk{8} = [ksai_trig(5),ksai_state(8)];
        
        Tk{i} = sort(Tk{i});
        R = [];
        
        % 筛选有效信息
        for j = 1 : 1 : length(Nk{i})
            if Nk{i}(j) >= Tk{i}(1) && Nk{i}(j) <= Tk{i}(end)
                R(end + 1) = Nk{i}(j);
            end
        end
        sigma = 0;
        for j = 1 : 1 : length(R)
            sigma = sigma + R(j);
        end
        sigma = sigma / length(R);
        if i <= 3
            ksai_state(i) = sigma + (P_axis{i}(end -1) -  ploss(end - 1)) - (P_axis{i}(end) -  ploss(end));
        else
            ksai_state(i) = sigma  + P_axis{i}(end) - P_axis{i}(end - 1);
        end
        
        % 对 lambda 进行攻击模拟
        ksai_state(6) = 200;
        ksai_state(7) = 400 * sin(0.004 * pi * k);
        ksai_state(8) = (k/30)^2 - 400;
        
        Ksai{i}(end + 1) = ksai_state(i);
    end
    
    %计算发电端与负载端的总功率差值
    error_temp = 0;
    for i = 1 : 1 : 3
        error_temp = error_temp + P_state(i);
    end
    error_temp = error_temp - pls;
    for i = 4 : 1 : 8
        error_temp = error_temp - P_state(i);
    end
    delta_P(end + 1) = error_temp;
    
    %事件触发判定
    for i = 1 : 1 : n
        error1 = abs(lambda_state(i) - lambda_trig(i));
        error2 = abs(ksai_state(i) - ksai_trig(i));
        if error1 >= Ek || error2 >= Ek
            lambda_trig(i) = lambda_state(i);
            ksai_trig(i) = ksai_state(i);
            count(i) = count(i) + 1;
            k_trig{i}(end + 1) = k;
        end
    end
    
    k_axis(end + 1) = k;  
end
end