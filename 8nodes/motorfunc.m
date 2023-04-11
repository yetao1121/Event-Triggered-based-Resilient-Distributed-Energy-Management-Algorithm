function P_star = motorfunc(i, lambda, p_m, p_M)

%优化发动机节点功率

a = [0.0024, 0.0056, 0.0072];
b = [5.56, 4.32, 6.60];
c = [30, 25, 25];
 
 
 f = @(x) a(i) * x^2 + b(i) * x + c(i) - lambda * x;
 
 P_star = fminbnd(f, p_m, p_M);
 
end