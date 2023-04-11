    function W = doubly_stochastic()

% 此双随机矩阵根据 IEEE 39-bus 系统拓扑结构计算来的权重矩阵，值得注意的是，我们考虑的是无向图

W = [0.642857142857143,0,0,0,0,0,0,0,0,0,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0,0,0,0,0;
    0,0.107142857142858,0.178571428571429,0,0,0,0,0,0,0,0,0,0.178571428571429,0.178571428571429,0.178571428571429,0,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0;
    0,0.178571428571429,0.642857142857143,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0.178571428571429,0.642857142857143,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0.178571428571429,0.464285714285715,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0.178571428571429,0.464285714285715,0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0.178571428571429,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0.464285714285715,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0.178571428571429,0,0;
    0,0,0,0,0,0,0.178571428571429,0.642857142857143,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0,0;
    0,0,0,0,0,0,0,0,0.464285714285715,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0,0,0,0.178571428571429;
    0,0,0,0,0,0,0,0,0.178571428571429,0.642857142857143,0,0,0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0,0,0,0;
    0.178571428571429,0,0,0,0,0,0,0,0,0,0.464285714285715,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0.178571428571429,0.464285714285715,0,0,0,0.178571428571429,0,0,0,0,0,0,0.178571428571429,0,0,0,0,0;
    0,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0.642857142857143,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0.642857142857143,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0,0.642857142857143,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0,0,0.178571428571429,0.464285714285715,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0;
    0,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0.464285714285715,0.178571428571429,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0.178571428571429,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0.464285714285715,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0,0.642857142857143,0,0.178571428571429,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.464285714285715,0.178571428571429,0,0,0,0.178571428571429,0,0.178571428571429,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0.178571428571429,0.642857142857143,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.464285714285715,0.178571428571429,0.178571428571429,0,0,0,0.178571428571429;
    0.178571428571429,0,0,0,0,0,0,0,0,0,0.178571428571429,0.178571428571429,0,0,0,0,0,0,0,0,0,0.178571428571429,0.285714285714286,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0.178571428571429,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0,0.285714285714286,0,0,0,0.178571428571429;
    0,0,0,0,0,0,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0,0,0,0,0.285714285714286,0.178571428571429,0.178571428571429,0;
    0,0,0,0,0,0,0.178571428571429,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0.464285714285715,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0,0,0,0,0.178571428571429,0,0.464285714285715,0.178571428571429;
    0,0,0,0,0,0,0,0,0.178571428571429,0,0,0,0,0,0,0,0,0,0,0,0,0.178571428571429,0,0.178571428571429,0,0,0.178571428571429,0.285714285714286];

end