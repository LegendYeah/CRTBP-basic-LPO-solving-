%% Differential correction (time_fix & time_free) in continunig Earth-Moon L2 halo orbit family 

clear
clc
Miu = 0.0121505482564457;%mass ratio
[X1, X2, X3, Gamma1, Gamma2, Gamma3] = LibrationPoints(Miu);
Gamma = Gamma2; 
X_i = X2;
opts = odeset('RelTol',1e-13,'AbsTol',1e-15);
Amplitude = 0.045; %  Amplitude
[X0, period] = HaloThirdOrder( Amplitude, 'Az', 'L2', Miu, 0, 0)
xx0_exact0 = X0';

 %% 01 DC_time_free method 


tic
[xx0_exact,T_DC,Omega_DC_timefree] = DC_timefree(xx0_exact0,Miu);
[~,XOrbit_DC_timefree] = ode113(@CR3BP, [0,1*T_DC],  xx0_exact, opts, Miu);
XOrbit_DC_timefree(end,:) - XOrbit_DC_timefree(1,:)



%% 02 DC_time_fix method

T_continue = T_DC;
initial_X_DC = xx0_exact;

for jj = 1:100

    initial_X_DC = DC_timefix(initial_X_DC, T_DC, Miu);
    [tOrbit_DC,XOrbit_DC] = ode113(@CR3BP,  [0,T_DC],  initial_X_DC , opts, Miu); 
    T_DC = T_DC-1e-3;
    figure(2)
    hold on
    y2 = plot3(XOrbit_DC(:,1),XOrbit_DC(:,2),XOrbit_DC(:,3),'LineWidth',1)
    
end

figure(2)
y1 = plot3(XOrbit_DC_timefree(:,1),XOrbit_DC_timefree(:,2),XOrbit_DC_timefree(:,3),'r','LineWidth',2) % initial seed orbit
hold on
% y2 = plot3(XOrbit_DC(:,1),XOrbit_DC(:,2),XOrbit_DC(:,3),'k','LineWidth',2)
grid on
box on
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
set(gca, 'FontSize',20,'Fontname', 'Times New Roman');






