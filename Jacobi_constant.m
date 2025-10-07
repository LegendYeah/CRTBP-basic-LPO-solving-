function C = Jacobi_constant(xx1,Miu)
%%%%%%%%% 1） 求解Jacobi constant, xx1 为质心会和坐标系下的归一化状态
    x = xx1(1); y = xx1(2); z = xx1(3);
    xdot = xx1(4); ydot = xx1(5); zdot = xx1(6);
    r1 = sqrt((x + Miu)^2 + y^2 + z^2);
    r2 = sqrt((x - 1 + Miu)^2 + y^2 + z^2);
    C = (x^2 + y^2 ) + 2 * (1 - Miu) / r1 + 2* Miu / r2 + Miu * (1 - Miu) - (xdot^2 + ydot^2 + zdot^2); % Jacobi 常数   
end
