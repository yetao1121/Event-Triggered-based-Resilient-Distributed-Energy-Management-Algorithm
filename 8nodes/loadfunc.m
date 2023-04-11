function P_star = loadfunc(i,k,lambda,p_m,p_M)

%优化负载节点功率

    w = [0, 0, 0, 18.43, 13.17, 15.46, 10.03, 8.45];
    e = [0, 0, 0, 0.0545, 0.0877, 0.0547, 0.1041, 0.0870];

    %f = @(x) lambda * x - (w(i) * x - e(i) * x^2);

    %g = @(x) lambda * x - ((w(i)^2) / (4 * e(i)));zz
    if k >= 150
        if i >= 6
            L = 6.26;
        else
            L = lambda;
        end
    else
        L = lambda;
    end
    if L < 0
        if(w(i)/(2 * e(i)) >= p_M)
            P_star = p_M;
        elseif (w(i)/(2 * e(i)) <= p_m)
            P_star = p_m;
        else
            P_star = w(i)/(2 * e(i));
        end
    else
        if (w(i) - L) / (2 * e(i)) > p_m
            if (w(i) - L) / (2 * e(i)) > p_M
                P_star = p_M;
            else
                P_star = (w(i) - L) / (2 * e(i));
            end
        else
            P_star = p_m;
        end
    end
end