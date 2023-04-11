[SW_CEMA, k_axis] = CEMA_SW();
[SW_ETRDEMA, ~] = ETRDEMA_SW();
plot(k_axis, -SW_CEMA, 'LineWidth', 2);
hold on;
grid on;
plot(k_axis, -SW_ETRDEMA, 'LineWidth', 2);
xlabel('k');
ylabel('社会福利');
legend('传统能量管理算法', '本章所提方法');