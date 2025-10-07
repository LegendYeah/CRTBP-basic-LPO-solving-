%% 三体周期轨道微分校正_ 积分时间固定

function [xx0_exact,counter] = DC_timefix(xx0, T, mu)  %% T 为积分时间，xx0 为6维状态量, 要求xx0 必须位于对xoz平面上，对称点
   
%     if isempty(Accuracy)
%         Accuracy = 5e-13;
%     end

    STM_I = reshape(eye(6),[],36);% Jacobi 矩阵初值
    xx0 = [xx0, STM_I];
    ode_options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    deltaVec = [100; 100; 100];
    counter=0;
    verbose = 1;

    while (abs(norm(deltaVec)) > 3e-15 && counter < 25)

        counter = counter+1;
        [tHalfOrbit,XHalfOrbit] = ode113(@CR3BP, [0 T/2], xx0, ode_options, mu);
    
        % State when y=0 (half orbit period)
        dX = XHalfOrbit(end,1:6);
        
        % STM at half orbit period
        dPhi = reshape(XHalfOrbit(end,7:end),6,[]); % 积分末端时刻的Jacobi 矩阵
    
        %% dX/dt (X is state)
        xdot = dX(4);
        ydot = dX(5);
        zdot = dX(6);
    
        r1 = sqrt((dX(1)+mu)^2 + dX(2)^2 + dX(3)^2); %S/C distance to Sun
        r2 = sqrt((dX(1)-1+mu)^2 + dX(2)^2 + dX(3)^2); %S/C distance to Earth
    
        % Accelerations
        xdotdot = 2*ydot + dX(1) - (1 - mu)*((dX(1) + mu)/(r1^3)) - mu*(dX(1) - 1 + mu)/(r2^3);
        ydotdot = -2*xdot + dX(2) - (1 - mu)*(dX(2)/(r1^3)) - mu*(dX(2))/(r2^3);
        zdotdot = -(1 - mu)*(dX(3))/(r1^3) - mu*(dX(3))/(r2^3); 
        
        % Derivative of new state
        dXdt = [xdot ydot zdot xdotdot ydotdot zdotdot];

        %解算x0和doty0,进行微分修正


        updateMat=[dPhi(2,1) dPhi(2,3) dPhi(2,5); dPhi(4,1) dPhi(4,3) dPhi(4,5); dPhi(6,1) dPhi(6,3) dPhi(6,5)];
        deltaVec = pinv(updateMat)*[-XHalfOrbit(end,2); -XHalfOrbit(end,4); -XHalfOrbit(end,6)];
        
        if verbose
            norm(deltaVec)
            fprintf('Iteration counter2: %d\n', counter);
        end

        if (abs(xx0(3)) <=1e-25 && abs(xx0(1) - (1-mu))<=1e-2) %%%%%%%%%%%%%%%%  针对halo轨道求解
            break;
        end

        deltaX = [deltaVec(1), 0, deltaVec(2), 0, deltaVec(3), 0];
    
        xx0(1:6) = xx0(1:6) + deltaX;

    end
     xx0_exact = xx0(1:6);

end