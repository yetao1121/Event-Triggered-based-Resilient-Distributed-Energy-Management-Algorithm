clear; clc;

n = 28;% 节点数量

% 发电端参数
a = [0.0024, 0.0056, 0.0072, 0.0047, 0.0091, 0.0018, 0.0053, 0.0063, 0.0028, 0.0046];
b = [5.56, 4.32, 6.60, 3.14, 7.54, 3.28, 7.31, 2.45, 7.63, 4.76];
c = [30, 25, 25, 16, 6, 54, 23, 15, 20, 12];

% 负载端参数
w = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18.43, 13.17, 15.46, 10.03, 8.45, 15.38, 19.16, 16.85, 15.63, 6.75, 14.95, 5.87, 18.18, 15.08, 14.90, 19.45, 18.05, 15.83];
e = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0545, 0.0877, 0.0547, 0.1041, 0.0870, 0.0984, 0.1564, 0.0564, 0.0950, 0.0470, 0.0970, 0.0349, 0.0879, 0.0653, 0.0897, 0.1345, 0.0924, 0.1026];

% 局部约束上下限
p_min = [60, 25, 28, 40, 35, 29, 45, 56, 12, 30, 50, 100, 40, 30, 80, 40, 80, 50, 50, 78, 103, 33, 99, 89, 123, 5, 43, 67];
p_max = [339.69, 479.10, 290.4, 306.31, 593.80, 137.19, 595.40, 162.17, 165.1, 443.41, 100.34, 159.13, 80.56, 123.98, 109.55, 76.34, 137.93, 84.19, 104.06, 119.36, 176.19, 88.65, 175.31, 129.03, 167.42, 326.51, 79.38, 147.26];

 %初始状态
 for i = 1 : 1 : n
     if i <= 10
         lambda_state(i) = (2 * a(i) * p_max(i) + b(i));
         lambda_trig(i) = lambda_state(i);
     else
         lambda_state(i) = -2 * e(i) * p_min(i) + w(i);
         lambda_trig(i) = lambda_state(i);
     end
 end
P_state = zeros(1,n);    % 初始化每个节点的功率
ksai_state = zeros(1,n); % 初始化每个节点的 ksai 值
ksai_trig = zeros(1, n);
k_max = 1000;
count = zeros(1,n);      % 记录每个节点的触发次数
count_notrig = k_max * ones(1, n); % 800次迭代，没有事件触发每个节点都是800次触发，做柱状图对比
eta = 0.000708;

%定义数组存放每一次迭代的数据以用于作图
for i = 1:1:n
    P_axis{i} = 0; 
    Lambda{i} = 0;
    ksai{i} = 0;
    k_trig{i} = []; 
end
k_axis  = 0;
ploss = 0;
delta_p = 0;
SW = 0;
optimal = 8300;

for k = 1 : 1 : k_max
    ak = 1 / (1 + 20 * k);   %衰减迭代步长
    Ek = 5 / (1 + 0.5 * k);  %事件触发阈值
    
    for i = 1 : 1 : n
        %定义包含自身以及邻居的状态集（其中邻居为上一次触发时刻状态，自身为实时状态）
        Nl{1} = [lambda_state(1), lambda_trig(10), lambda_trig(11), lambda_trig(23)];
        Nl{2} = [lambda_state(2), lambda_trig(3), lambda_trig(13), lambda_trig(14), lambda_trig(15), lambda_trig(17)];
        Nl{3} = [lambda_trig(2), lambda_state(3), lambda_trig(4)];
        Nl{4} = [lambda_trig(3), lambda_state(4), lambda_trig(5), lambda_trig(17)];
        Nl{5} = [lambda_trig(4), lambda_state(5), lambda_trig(6), lambda_trig(18)];
        Nl{6} = [lambda_trig(5), lambda_state(6), lambda_trig(18), lambda_trig(19)];
        Nl{7} = [lambda_state(7), lambda_trig(8), lambda_trig(19), lambda_trig(25), lambda_trig(26)];
        Nl{8} = [lambda_trig(7), lambda_state(8), lambda_trig(26)];
        Nl{9} = [lambda_state(9), lambda_trig(10), lambda_trig(24), lambda_trig(28)];
        Nl{10} = [lambda_trig(1), lambda_trig(9), lambda_state(10), lambda_trig(24)];
        Nl{11} = [lambda_trig(1), lambda_state(11), lambda_trig(12), lambda_trig(23)];
        Nl{12} = [lambda_trig(11), lambda_state(12), lambda_trig(14), lambda_trig(16), lambda_trig(21), lambda_trig(23)];
        Nl{13} = [lambda_trig(2), lambda_state(13), lambda_trig(14)];
        Nl{14} = [lambda_trig(2), lambda_trig(12), lambda_trig(13), lambda_state(14)];
        Nl{15} = [lambda_trig(2), lambda_state(15), lambda_trig(16)];
        Nl{16} = [lambda_trig(12), lambda_trig(15), lambda_state(16), lambda_trig(17), lambda_trig(21)];
        Nl{17} = [lambda_trig(2), lambda_trig(4), lambda_trig(16), lambda_state(17), lambda_trig(18)];
        Nl{18} = [lambda_trig(5), lambda_trig(6), lambda_trig(17), lambda_state(18), lambda_trig(19)];
        Nl{19} = [lambda_trig(6), lambda_trig(7), lambda_trig(18), lambda_state(19), lambda_trig(21)];
        Nl{20} = [lambda_state(20), lambda_trig(21), lambda_trig(25), lambda_trig(27)];
        Nl{21} = [lambda_trig(12), lambda_trig(16), lambda_trig(19), lambda_trig(20), lambda_state(21)];
        Nl{22} = [lambda_state(22), lambda_trig(23), lambda_trig(24), lambda_trig(28)];
        Nl{23} = [lambda_trig(1), lambda_trig(11), lambda_trig(12), lambda_trig(22), lambda_state(23)];
        Nl{24} = [lambda_trig(9), lambda_trig(10), lambda_trig(22), lambda_state(24), lambda_trig(28)];
        Nl{25} = [lambda_trig(7), lambda_trig(20), lambda_state(25), lambda_trig(26), lambda_trig(27)];
        Nl{26} = [lambda_trig(7), lambda_trig(8), lambda_trig(25), lambda_state(26), lambda_trig(27)];
        Nl{27} = [lambda_trig(20), lambda_trig(25), lambda_trig(26), lambda_state(27), lambda_trig(28)];
        Nl{28} = [lambda_trig(9), lambda_trig(22), lambda_trig(24), lambda_trig(27), lambda_state(28)];

        %定义包含自身以及可信邻居的状态集（其中邻居为上一次触发时刻状态，自身为实时状态）
        Tl{1} = [lambda_state(1), lambda_trig(23)];
        Tl{2} = [lambda_state(2), lambda_trig(17)];
        Tl{3} = [lambda_trig(2), lambda_state(3)];
        Tl{4} = [lambda_state(4), lambda_trig(17)];
        Tl{5} = [lambda_state(5), lambda_trig(18)];
        Tl{6} = [lambda_state(6), lambda_trig(18), lambda_trig(19)];
        Tl{7} = [lambda_state(7), lambda_trig(19)];
        Tl{8} = [lambda_trig(7), lambda_state(8)];
        Tl{9} = [lambda_state(9), lambda_trig(28)];
        Tl{10} = [lambda_trig(1), lambda_state(10)];
        Tl{11} = [lambda_trig(1), lambda_state(11), lambda_trig(12), lambda_trig(23)];
        Tl{12} = [lambda_state(12),lambda_trig(21), lambda_trig(23)];
        Tl{13} = [lambda_trig(2), lambda_state(13)];
        Tl{14} = [lambda_trig(2), lambda_trig(12), lambda_state(14)];
        Tl{15} = [lambda_trig(2), lambda_state(15)];
        Tl{16} = [lambda_trig(12), lambda_state(16), lambda_trig(17), lambda_trig(21)];
        Tl{17} = [lambda_trig(2), lambda_state(17)];
        Tl{18} = [lambda_trig(17), lambda_state(18), lambda_trig(19)];
        Tl{19} = [lambda_trig(7), lambda_trig(18), lambda_state(19), lambda_trig(21)];
        Tl{20} = [lambda_state(20), lambda_trig(21)];
        Tl{21} = [lambda_trig(12), lambda_trig(19), lambda_state(21)];
        Tl{22} = [lambda_state(22), lambda_trig(23), lambda_trig(28)];
        Tl{23} = [lambda_trig(1), lambda_trig(12), lambda_trig(22), lambda_state(23)];
        Tl{24} = [lambda_trig(22), lambda_state(24), lambda_trig(28)];
        Tl{25} = [lambda_trig(7), lambda_state(25)];
        Tl{26} = [lambda_trig(7), lambda_state(26)];
        Tl{27} = [lambda_state(27), lambda_trig(28)];
        Tl{28} = [lambda_trig(22), lambda_state(28)];
        
        Tl{i} = sort(Tl{i});
        R = [];
        % 进行筛选
        for j = 1 : 1 : length(Nl{i})
            if Nl{i}(j) <= Tl{i}(end) && Nl{i}(j) >= Tl{i}(1)
                R(end + 1) = Nl{i}(j);
            end
        end
        
        sigma = 0;
        for j = 1 : 1 : length(R)
            sigma = sigma + R(j);
        end
        sigma = sigma / length(R);
        lambda_state(i) = sigma + eta * ksai_state(i);
        
        % 对 lambda 进行攻击模拟
        lambda_state(11) = 10;
        lambda_state(13) = 8 * sin(0.003 * pi * k) + 6;
        lambda_state(15) = (k/300)^2;
        lambda_state(25) = -(k/300)^2 + 12;
        lambda_state(26) = -1;
        
        Lambda{i}(end+1) = lambda_state(i);
    end
    
    % 更新功率
    for i = 1 : 1 : n
        if i <= 10
             P_state(i) = motorfunc(i, lambda_state(i), p_min(i), p_max(i));
        else
            P_state(i) = loadfunc(i, k, lambda_state(i), p_min(i), p_max(i));
        end
        P_axis{i}(end+1) = P_state(i);
    end
    
    % 计算目标函数值
    sw = 0;
    for i = 1 : 1 : 10
        sw = sw + a(i) * P_state(i)^2 + b(i) * P_state(i) + c(i);
    end
    for i = 11 : 1 : 28
        if P_state(i) <= w(i)/(2 * e(i))
            sw = sw - w(i) * P_state(i) + e(i) * P_state(i);
        else
            sw = sw - w(i)^2/(4 * e(i));
        end
    end
    SW(end + 1) = sw + 3200;
    optimal(end + 1) = 8300;
    
    %更新 ksai
    for i = 1 : 1 : n
        %定义包含自身以及邻居的状态集（其中邻居为上一次触发时刻状态，自身为实时状态）
        Nk{1} = [ksai_state(1), ksai_trig(10), ksai_trig(11), ksai_trig(23)];
        Nk{2} = [ksai_state(2), ksai_trig(3), ksai_trig(13), ksai_trig(14), ksai_trig(15), ksai_trig(17)];
        Nk{3} = [ksai_trig(2), ksai_state(3), ksai_trig(4)];
        Nk{4} = [ksai_trig(3), ksai_state(4), ksai_trig(5), ksai_trig(17)];
        Nk{5} = [ksai_trig(4), ksai_state(5), ksai_trig(6), ksai_trig(18)];
        Nk{6} = [ksai_trig(5), ksai_state(6), ksai_trig(18), ksai_trig(19)];
        Nk{7} = [ksai_state(7), ksai_trig(8), ksai_trig(19), ksai_trig(25), ksai_trig(26)];
        Nk{8} = [ksai_trig(7), ksai_state(8), ksai_trig(26)];
        Nk{9} = [ksai_state(9), ksai_trig(10), ksai_trig(24), ksai_trig(28)];
        Nk{10} = [ksai_trig(1), ksai_trig(9), ksai_state(10), ksai_trig(24)];
        Nk{11} = [ksai_trig(1), ksai_state(11), ksai_trig(12), ksai_trig(23)];
        Nk{12} = [ksai_trig(11), ksai_state(12), ksai_trig(14), ksai_trig(16), ksai_trig(21), ksai_trig(23)];
        Nk{13} = [ksai_trig(2), ksai_state(13), ksai_trig(14)];
        Nk{14} = [ksai_trig(2), ksai_trig(12), ksai_trig(13), ksai_state(14)];
        Nk{15} = [ksai_trig(2), ksai_state(15), ksai_trig(16)];
        Nk{16} = [ksai_trig(12), ksai_trig(15), ksai_state(16), ksai_trig(17), ksai_trig(21)];
        Nk{17} = [ksai_trig(2), ksai_trig(4), ksai_trig(16), ksai_state(17), ksai_trig(18)];
        Nk{18} = [ksai_trig(5), ksai_trig(6), ksai_trig(17), ksai_state(18), ksai_trig(19)];
        Nk{19} = [ksai_trig(6), ksai_trig(7), ksai_trig(18), ksai_state(19), ksai_trig(21)];
        Nk{20} = [ksai_state(20), ksai_trig(21), ksai_trig(25), ksai_trig(27)];
        Nk{21} = [ksai_trig(12), ksai_trig(16), ksai_trig(19), ksai_trig(20), ksai_state(21)];
        Nk{22} = [ksai_state(22), ksai_trig(23), ksai_trig(24), ksai_trig(28)];
        Nk{23} = [ksai_trig(1), ksai_trig(11), ksai_trig(12), ksai_trig(22), ksai_state(23)];
        Nk{24} = [ksai_trig(9), ksai_trig(10), ksai_trig(22), ksai_state(24), ksai_trig(28)];
        Nk{25} = [ksai_trig(7), ksai_trig(20), ksai_state(25), ksai_trig(26), ksai_trig(27)];
        Nk{26} = [ksai_trig(7), ksai_trig(8), ksai_trig(25), ksai_state(26), ksai_trig(27)];
        Nk{27} = [ksai_trig(20), ksai_trig(25), ksai_trig(26), ksai_state(27), ksai_trig(28)];
        Nk{28} = [ksai_trig(9), ksai_trig(22), ksai_trig(24), ksai_trig(27), ksai_state(28)];

        %定义包含自身以及可信邻居的状态集（其中邻居为上一次触发时刻状态，自身为实时状态）
        Tk{1} = [ksai_state(1), ksai_trig(23)];
        Tk{2} = [ksai_state(2), ksai_trig(17)];
        Tk{3} = [ksai_trig(2), ksai_state(3)];
        Tk{4} = [ksai_state(4), ksai_trig(17)];
        Tk{5} = [ksai_state(5), ksai_trig(18)];
        Tk{6} = [ksai_state(6), ksai_trig(18), ksai_trig(19)];
        Tk{7} = [ksai_state(7), ksai_trig(19)];
        Tk{8} = [ksai_trig(7), ksai_state(8)];
        Tk{9} = [ksai_state(9), ksai_trig(28)];
        Tk{10} = [ksai_trig(1), ksai_state(10)];
        Tk{11} = [ksai_trig(1), ksai_state(11), ksai_trig(12), ksai_trig(23)];
        Tk{12} = [ksai_state(12),ksai_trig(21), ksai_trig(23)];
        Tk{13} = [ksai_trig(2), ksai_state(13)];
        Tk{14} = [ksai_trig(2), ksai_trig(12), ksai_state(14)];
        Tk{15} = [ksai_trig(2), ksai_state(15)];
        Tk{16} = [ksai_trig(12), ksai_state(16), ksai_trig(17), ksai_trig(21)];
        Tk{17} = [ksai_trig(2), ksai_state(17)];
        Tk{18} = [ksai_trig(17), ksai_state(18), ksai_trig(19)];
        Tk{19} = [ksai_trig(7), ksai_trig(18), ksai_state(19), ksai_trig(21)];
        Tk{20} = [ksai_state(20), ksai_trig(21)];
        Tk{21} = [ksai_trig(12), ksai_trig(19), ksai_state(21)];
        Tk{22} = [ksai_state(22), ksai_trig(23), ksai_trig(28)];
        Tk{23} = [ksai_trig(1), ksai_trig(12), ksai_trig(22), ksai_state(23)];
        Tk{24} = [ksai_trig(22), ksai_state(24), ksai_trig(28)];
        Tk{25} = [ksai_trig(7), ksai_state(25)];
        Tk{26} = [ksai_trig(7), ksai_state(26)];
        Tk{27} = [ksai_state(27), ksai_trig(28)];
        Tk{28} = [ksai_trig(22), ksai_state(28)];
        
        Tk{i} = sort(Tk{i});
        R = [];
        % 进行筛选
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
        if i <= 10
            ksai_state(i) = sigma + (P_axis{i}(end -1)) - (P_axis{i}(end));
        else
            ksai_state(i) = sigma  + P_axis{i}(end) - P_axis{i}(end - 1);
        end
        % 对 ksai 进行攻击模拟
        ksai_state(11) = 50;
        ksai_state(13) = 80 * sin(0.003 * pi * k);
        ksai_state(15) = (k/100)^2-20;
        ksai_state(25) = -(k/100)^2 + 60;
        ksai_state(26) = -50;
        
        ksai{i}(end+1) = ksai_state(i);
    end
    
    %计算发电端与负载端的总功率差值
    error_temp = 0;
    for i = 1 : 1 : 10
        error_temp = error_temp + P_state(i);
    end
    error_temp = error_temp;
    for j = 11 : 1 : 28
        error_temp = error_temp - P_state(j);
    end
    delta_p(end + 1) = error_temp + 500;
    
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

% figure(1); % lambda的迭代过程
subplot(2, 2, 1);
for i = 1 : 1 : n
    if i == 11 || i == 13 || i == 15 || i == 25 || i == 26
        plot(k_axis,Lambda{i},'--','LineWidth',1.5);
    else
        plot(k_axis,Lambda{i},'LineWidth',1.5);
    end
    hold on;
end
grid on;
xlabel('迭代次数 k');
ylabel('\lambda_i(k)');

% P的迭代过程
subplot(2, 2, 2);
for i = 1:1:10
    plot(k_axis,P_axis{i},'LineWidth',1.5);
    hold on;
end
for i = 11:1:28
    plot(k_axis,-P_axis{i},'LineWidth',1.5);
    hold on;
    grid on;
end
xlabel('迭代次数 k');
ylabel('输出功率 P_i(MW)');

% ksai的迭代过程
subplot(2, 2, 3);
for i = 1:1:n
    if i == 11 || i == 13 || i == 15 || i == 25 || i == 26
        plot(k_axis,ksai{i},'--','LineWidth',1.5);
    else
        plot(k_axis,ksai{i},'LineWidth',1.5);
    end
    hold on;
end
grid on;
axis([0 1000 -100 100]);
xlabel('迭代次数 k');
ylabel('偏差电量 \xi_i');

% 负载平衡约束的动态过程
subplot(2, 2, 4);
plot(k_axis, -delta_p, 'LineWidth', 2);
xlabel('迭代次数 k');
ylabel('总功率平衡约束(MW)');
grid on;

figure(2); % 事件触发通信次数对比
subplot(2, 1, 1);
x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28];
y = [count',count_notrig'];
bar(x,y);
hold on;
xlabel('节点编号');
ylabel('触发次数');
legend('事件触发通信','传统周期通信');

subplot(2, 1, 2);
tmp_node = 2;
trig_value = [];
for i = 1:1:length(k_trig{tmp_node})
    trig_value(end + 1) = Lambda{tmp_node}(k_trig{tmp_node}(i));
end
stem(k_trig{tmp_node},trig_value,'LineWidth',0.5);
axis([0,1000,0,8]);
xlabel('迭代次数 k');
ylabel('触发时刻状态');

figure(3)
SW_CEMA = CEMA_28nodes();
plot(k_axis, -SW, 'LineWidth', 2);
hold on;
plot(k_axis, -SW_CEMA, 'LineWidth', 2);
plot(k_axis, optimal, 'LineWidth', 2);
xlabel('k');
ylabel('社会福利');
legend('本章所提方法','未遭受攻击的能量管理方法','最优值');
grid on; 