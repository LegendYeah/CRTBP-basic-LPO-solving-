function [xx1,T,omega,counter] = DC_timefree(xx0,mu)

        ode_options = odeset('Events',@Findyzero,'RelTol',1e-13,'AbsTol',1e-15);%积分终止条件
        deltaVec = [100; 100];
        counter=0;
        verbose = 1;
        STM_I = reshape(eye(6),[],36);% Jacobi 矩阵初值
        xx0 = [xx0,STM_I] ;

        % 微分校正过程

        while ( abs(norm(deltaVec)) > 5e-14 && counter<=30)
            counter = counter+1;
            [tHalforbit,XHalforbit] = ode45(@CR3BP, [0 Inf], xx0, ode_options, mu);
        
            % State when y=0 (half orbit period)
            dX = XHalforbit(end,1:6);
            
            % STM at half orbit period
            dPhi = reshape(XHalforbit(end,7:end),6,[]); % 积分末端时刻的Jacobi 矩阵
        
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
            
            updateMat=[dPhi(2,1) dPhi(2,5) ydot; dPhi(4,1) dPhi(4,5) xdotdot; dPhi(6,1) dPhi(6,5) zdotdot];
            deltaVec = pinv(updateMat)*[0; -XHalforbit(end,4); -XHalforbit(end,6)];
        
        %     updateMat=[dPhi(2,3) dPhi(2,5) ydot; dPhi(4,3) dPhi(4,5) xdotdot; dPhi(6,3) dPhi(6,5) zdotdot];
        %     deltaVec = pinv(updateMat)*[0; -XHalforbit(end,4); -XHalforbit(end,6)];
            
            
            if verbose
                norm(deltaVec);
                fprintf('Iteration counter1: %d\n', counter)
            end
        
            deltaX = [deltaVec(1), 0, 0, 0, deltaVec(2), 0];
        %  deltaX = [0, 0, deltaVec(1), 0, deltaVec(2), 0];
            xx0(1:6) = xx0(1:6) + deltaX;

        end
        xx1 = xx0(1:6);
        T = 2*tHalforbit(end);
        omega = 2*pi/T; 
  