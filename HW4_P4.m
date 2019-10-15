% Gurpal Singh
% 10/27/2017
% Math 567 HW4

% Problem 4
clear all; close all; clc;

% Exact solution we derived in class
u = @(t) 1/4 * (2*t - 1 + 5*exp(-2*t));
f = @(t,u) t - 2*u; % RHS

% Convergence Study k values
kvec = [0.5 0.3 0.1 0.05 0.03 0.01 0.001 0.0005];
fprintf("********************************************************\n");
fprintf("\t \t \t \t \t  \t Error\n");
fprintf("********************************************************\n");
fprintf("k \t\t\t Exact \t\t\t\t FE \t\t\t\t RK4\n");
fprintf("--------------------------------------------------------\n");

for j = 1:length(kvec)
    % Time steps size
    k = kvec(j);
    
    % Final Time
    T = 2;
    
    % Vector of time increments
    t = 0:k:T;
    
    % Exact Solution
    exact = u(t); % Exact Solution
    
    % Starting values using exact solution
    y(1) = exact(1);
    y(2) = exact(2);
    y(3) = exact(3);
    y(4) = exact(4);
    
    % Loop through AB4 w/starting values
    for i = 5:length(t)
        y(i) = y(i-1) + k *( (55/24)*f(t(i-1),y(i-1)) - (59/24)*f(t(i-2),y(i-2)) + (37/24)*f(t(i-3),y(i-3)) - (3/8)*f(t(i-4),y(i-4)));
    end
    
    
    
    % Using Forward Euler for starting value
    % Initial Condition
    u_FE = zeros(length(t),1);
    u_FE(1) = 1;
    
    % Loop to get 4 starting values using Forward Euler
    for i = 2:4
        u_FE(i) = u_FE(i-1) + k*f(t(i-1),u_FE(i-1));
    end
    
    % Loop through AB4 w/starting values from FE
    for i = 5:length(t)
        u_FE(i) = u_FE(i-1) + k *( (55/24)*f(t(i-1),u_FE(i-1)) - (59/24)*f(t(i-2),u_FE(i-2)) + (37/24)*f(t(i-3),u_FE(i-3)) - (3/8)*f(t(i-4),u_FE(i-4)));
    end
    
   
    % Using the 4th Order Runge-Kutta Method for starting vals
    u_RK4 = zeros(length(t),1);
    u_RK4(1) = 1;
    
    for i=2:4
        k1 = f(t(i-1),u_RK4(i-1));
    
        k2 = f(t(i-1) + k/2, y(i-1) + (k/2)*k1);
    
        k3 = f(t(i-1) + k/2, y(i-1) + (k/2)*k2);
    
        k4 = f(t(i-1) + k, y(i-1) + k*k3);
        
        u_RK4(i) = u_RK4(i-1) + (k/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
    
    % Loop through AB4 w/starting values from FE
    for i = 5:length(t)
        u_RK4(i) = u_RK4(i-1) + k *( (55/24)*f(t(i-1),u_RK4(i-1)) - (59/24)*f(t(i-2),u_RK4(i-2)) + (37/24)*f(t(i-3),u_RK4(i-3)) - (3/8)*f(t(i-4),u_RK4(i-4)));
    end
    
    
    % Error
    y_error(j) = abs( exact(end) - y(end) );
    
    u_FE_error(j) = abs( exact(end) - u_FE(end));
    
    u_RK4_error(j) = abs( exact(end) - u_RK4(end));
    
    fprintf("%.4f \t %12.4e\t \t %12.4e \t %12.4e\n",k,y_error(j),u_FE_error(j),u_RK4_error(j));
    
    
    % Create plot
    if (k == 0.05)
        figure
        plot(t,exact,'black')
        hold on
        plot(t,y,'b*')
        hold on
        plot(t,u_FE,'r--')
        hold on
        plot(t,u_RK4,'go')
        hold on
        legend('Exact Solution','Exact Starting Values AB4','Forward Euler AB4','Runge-Kutta AB4')
        title('k = 0.05')
        xlabel('Time t')
        ylabel('Solution U(t)')
    end
    
end


for j = 2:length(kvec)
    Exact_Convergence_Rate(j-1) = (log(y_error(j)) - log(y_error(j-1)))  / (log(kvec(j)) - log(kvec(j-1)));
    
    FE_Convergence_Rate(j-1) = (log(u_FE_error(j)) - log(u_FE_error(j-1)))  / (log(kvec(j)) - log(kvec(j-1)));
    
    RK4_Convergence_Rate(j-1) = (log(u_RK4_error(j)) - log(u_RK4_error(j-1)))  / (log(kvec(j)) - log(kvec(j-1)));
end


% Printing Convergence Rate Results
% Only print last value 
fprintf("\n\n*************** The convergence rate for each of the methods ***************\n\n");
fprintf("Using Exact Starting Values: %.4f\n",Exact_Convergence_Rate(end));
fprintf("Using Forward Euler for Starting Values: %.4f\n",FE_Convergence_Rate(end));
fprintf("Using RK4 Starting Values: %.4f\n",RK4_Convergence_Rate(end));

Exact_Convergence_Rate
FE_Convergence_Rate
RK4_Convergence_Rate








