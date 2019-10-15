function [ u t ] = myFun( k, M, m )
%UNTITLED3 Takes in times step size k and Number of Timesteps M
%   Detailed explanation goes here

f = @(t,u) t - 2*u;    % RHS

T = 3;                 % Final time

uexact = @(t) 1/4 * (2*t - 1 + 5*exp(-2*t));        % Exact solution

u0 = uexact(0); 

% ------------------------------
% Solve the ODE
% ------------------------------

t = zeros(M+1,1);
u = zeros(M+1,1);

u(1) = u0;
t(1) = 0;   % = t0
method = m;       % 'fe', 'be', 'ab4', 'rk4'

for n = 1:M
    t(n+1) = n*k;
    switch method
        case 'exact'
            u(n+1) = uexact(t(n+1));
        case 'fe'
            u(n+1) = u(n) + k*f(t(n),u(n));
        case 'be'
            u(n+1) = 1/(1 + 2*k)*(u(n) + k*t(n+1));
        case 'heun'
            Y1 = u(n); 
            Y2 = u(n) + k*f(t(n),Y1);
            u(n+1) = u(n) + k/2*(f(t(n),Y1) + f(t(n) + k,Y2));
        case 'rk4'
            % 4th order Runge-Kutta method here
        case 'ab4'
            % Adam's Bashforth                
    end
end






end

