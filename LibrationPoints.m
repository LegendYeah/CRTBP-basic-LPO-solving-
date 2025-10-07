function [X2, X1, X3, Gamma1, Gamma2, Gamma3] = LibrationPoints(Miu)
    h = 0.01;
    error0 = 10;
    errorlim = 1e-10;
    X1 = 0;
    X2 = 0;
    X3 = 0;
    
    error = error0;
    while error > errorlim
        X1_Former = X1;
        K1 = -X1 + (1 - Miu) / (X1 + Miu) ^ 2 + Miu / (1 - X1 - Miu) ^ 2;
        K2 = -(X1 + h / 2 * K1) + (1 - Miu) / ((X1 + h / 2 * K1) + Miu) ^ 2 + Miu / (1 - (X1 + h / 2 * K1) - Miu) ^ 2;
        K3 = -(X1 + h / 2 * K2) + (1 - Miu) / ((X1 + h / 2 * K2) + Miu) ^ 2 + Miu / (1 - (X1 + h / 2 * K2) - Miu) ^ 2;
        K4 = -(X1 + h * K3) + (1 - Miu) / ((X1 + h * K3) + Miu) ^ 2 + Miu / (1 - (X1 + h * K3) - Miu) ^ 2;
        X1 = X1_Former + h / 6 * (K1 + 2 * K2 + 2 * K3+ K4);
        error = abs(X1 - X1_Former);
    end
    
    error = error0;
    while error > errorlim
        X2_Former = X2;
        K1 = -X2 + (1 - Miu) / (X2 + Miu) ^ 2 - Miu / (1 - X2 - Miu) ^ 2;
        K2 = -(X2 + h / 2 * K1) + (1 - Miu) / ((X2 + h / 2 * K1) + Miu) ^ 2 - Miu / (1 - (X2 + h / 2 * K1) - Miu) ^ 2;
        K3 = -(X2 + h / 2 * K2) + (1 - Miu) / ((X2 + h / 2 * K2) + Miu) ^ 2 - Miu / (1 - (X2 + h / 2 * K2) - Miu) ^ 2;
        K4 = -(X2 + h * K3) + (1 - Miu) / ((X2 + h * K3) + Miu) ^ 2 - Miu / (1 - (X2 + h * K3) - Miu) ^ 2;
        X2 = X2_Former + h / 6 * (K1 + 2 * K2 + 2 * K3+ K4);
        error = abs(X2 - X2_Former);
    end
    
    error = error0;
    while error > errorlim
        X3_Former = X3;
        K1 = -X3 - (1 - Miu) / (X3 + Miu) ^ 2 - Miu / (1 - X3 - Miu) ^ 2;
        K2 = -(X3 + h / 2 * K1) - (1 - Miu) / ((X3 + h / 2 * K1) + Miu) ^ 2 - Miu / (1 - (X3 + h / 2 * K1) - Miu) ^ 2;
        K3 = -(X3 + h / 2 * K2) - (1 - Miu) / ((X3 + h / 2 * K2) + Miu) ^ 2 - Miu / (1 - (X3 + h / 2 * K2) - Miu) ^ 2;
        K4 = -(X3 + h * K3) - (1 - Miu) / ((X3 + h * K3) + Miu) ^ 2 - Miu / (1 - (X3 + h * K3) - Miu) ^ 2;
        X3 = X3_Former + h / 6 * (K1 + 2 * K2 + 2 * K3+ K4);
        error = abs(X3 - X3_Former);
    end
    
    Gamma1 = 1 - Miu - X2;
    Gamma2 = X1 - 1 + Miu;
    Gamma3 = -Miu - X3;
end