function SW = CEMA_28nodes()

n = 28;% 节点数量
W = doubly_stochastic();
Q = W';

a = [0.0024, 0.0056, 0.0072, 0.0047, 0.0091, 0.0018, 0.0053, 0.0063, 0.0028, 0.0046];
b = [5.56, 4.32, 6.60, 3.14, 7.54, 3.28, 7.31, 2.45, 7.63, 4.76];
c = [30, 25, 25, 16, 6, 54, 23, 15, 20, 12];
B = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
w = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18.43, 13.17, 15.46, 10.03, 8.45, 15.38, 19.16, 16.85, 15.63, 6.75, 14.95, 5.87, 18.18, 15.08, 14.90, 19.45, 18.05, 15.83];
e = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0545, 0.0877, 0.0547, 0.1041, 0.0870, 0.0984, 0.1564, 0.0564, 0.0950, 0.0470, 0.0970, 0.0349, 0.0879, 0.0653, 0.0897, 0.1345, 0.0924, 0.1026];
p_min = [60, 25, 28, 40, 35, 29, 45, 56, 12, 30, 50, 100, 40, 30, 80, 40, 80, 50, 50, 78, 103, 33, 99, 89, 123, 5, 43, 67];
p_max = [339.69, 479.10, 290.4, 306.31, 593.80, 137.19, 595.40, 162.17, 165.1, 443.41, 100.34, 159.13, 80.56, 123.98, 109.55, 76.34, 137.93, 84.19, 104.06, 119.36, 176.19, 88.65, 175.31, 129.03, 167.42, 326.51, 79.38, 147.26];

  
%初始状态
lambda_state = [];
for i = 1:1:28
    if i<= 10
        lambda_state(end + 1) = (2 * a(i) * p_min(i) + b(i)) / (1 - 2 * B(i) * p_min(i));
    else
        if (p_max(i) > (w(i) / (2 * e(i))))
            lambda_state(end + 1) = 0;
        else
            lambda_state(end + 1) = w(i) - 2 * e(i) * p_max(i);
        end
    end
end

P_state = zeros(1,n); % 初始化每个节点的功率
ksai_state = zeros(1,n); 
eta = 0.00151;
k_max =1000;
SW = 0;

%定义数组存放每一次迭代的数据以用于作图
for i = 1:1:28
    P_axis{i} = 0; 
    Lambda{i} = 0;
    ksai{i} = 0;
end
k_axis  = 0;
ploss = 0;
error = 0; %#ok<*NBRAK>

for k = 1:1:k_max
    
    % 节点进行 lambda 迭代
    for i = 1:1:28
        sigma = 0;
        for j = 1:1:28
            sigma = sigma + W(i,j) * lambda_state(j);
        end
        lambda_state(i) = sigma + eta * ksai_state(i); 
        Lambda{i}(end+1) = lambda_state(i);
    end
    
    % 更新功率
    for i = 1:1:28
        if i <= 10
            P_state(i) = motorfunc_cema(i,lambda_state(i),p_min(i),p_max(i));
        else
            P_state(i) = loadfunc_cema(i,lambda_state(i),p_min(i),p_max(i));
        end
        P_axis{i}(end+1) = P_state(i);
    end
    
    % 计算目标函数值
    sw = 0;
    for i = 1:1:10
        sw = sw + a(i) * P_state(i)^2 + b(i) * P_state(i) + c(i);
    end
    for i = 11:1:28
        if P_state(i) <= w(i)/(2 * e(i))
            sw = sw - w(i) * P_state(i) + e(i) * P_state(i);
        else
            sw = sw - w(i)^2/(4 * e(i));
        end
    end
    SW(end + 1) = sw;
    
    %更新 ksai
    for i = 1:1:28
        sigma1 = 0;
        for j = 1:1:28
            sigma1 = sigma1 + Q(i,j) * ksai_state(j);
        end
        if i <= 10
            ksai_state(i) = sigma1 + (P_axis{i}(end -1) -  B(i) * (P_axis{i}(end -1))^2) - (P_axis{i}(end) -  B(i) * (P_axis{i}(end))^2);
        else
            ksai_state(i) = sigma1  + P_axis{i}(end) - P_axis{i}(end - 1);
        end
        ksai{i}(end+1) = ksai_state(i);
    end
    
    %计算发电端与负载端的总功率差值
    error_temp = 0;
    for i = 1:10
        error_temp = error_temp + P_state(i) - B(i) * P_state(i)^2;
    end
    for j = 11:28
        error_temp = error_temp - P_state(j);
    end
    error(end + 1) = error_temp;
    
    k_axis(end + 1) = k; 
end
end