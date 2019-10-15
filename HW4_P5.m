% Gurpal Singh
% 11/9/2017
% Math 567 HW4

% Problem 5

% Stability Region for AB4
rho_AB4 = @(zeta) zeta.^4 - zeta.^3;
sigma_AB4 = @(zeta)(-9/24) + (37/24)*zeta - (59/24)*zeta.^2 + (55/24)*zeta.^3;
z_AB4 = @(rho_AB4,sigma_AB4) rho_AB4./sigma_AB4;

theta = 0:.01:2*pi;
zeta = exp(i*theta);

rho_out = rho_AB4(zeta);
sigma_out = sigma_AB4(zeta);

z_out = z_AB4(rho_out,sigma_out);
ABz = z_out;

plot(z_out)
fill(real(z_out),imag(z_out),'r')
hold on
title('Stability Region for 4th order Adams Bashforth')
xlabel('Re(z)')
ylabel('Im(z)')
grid
axis([-1 1 -1 1])


% Stability Region for AM4
clear rho_out sigma_out z_out
rho_AM4 = @(zeta) zeta.^4 - zeta.^3;
sigma_AM4 = @(zeta) (1/720) * (-19 + 106*zeta - 264*zeta.^2 + 646*zeta.^3 + 251*zeta.^4);
z_AM4 = @(rho_AM4,sigma_AM4) rho_AM4./sigma_AM4;

rho_out = rho_AM4(zeta);
sigma_out = sigma_AM4(zeta);
z_out = z_AM4(rho_out,sigma_out);
AMz = z_out;
figure
plot(z_out)
hold on
fill(real(z_out),imag(z_out),'r')
title('Stability Region for 4th order Adams Moulton')
xlabel('Re(z)')
ylabel('Im(z)')
grid
axis([-3 1 -3 3])

% Showing RK4 Stability Region
figure
hold on
axis([-5 5 -5 5])
grid 

a = zeros(4,length(theta));

for k = 1:length(theta)
    c = [1./24 1./6 0.5 1 1-exp(i*theta(k))];
    a(:,k) = roots(c);
end

plot(a(1,:),'co:')
plot(a(2,:),'co:')
plot(a(3,:),'co:')
plot(a(4,:),'co:')
title('Stability Region for 4th order Runge-Kutta')
xlabel('Re(z)')
ylabel('Im(z)')
hold off

% Overlaying Plots to Compare
% AB4 + RK4

% RK4
figure
plot(a(1,:),'co:')
hold on
plot(a(2,:),'co:')
hold on
plot(a(3,:),'co:')
hold on
plot(a(4,:),'co:')
hold on


plot(ABz)
fill(real(ABz),imag(ABz),'r')
hold on
title('Stability Region for AB4 + RK4')
xlabel('Re(z)')
ylabel('Im(z)')
grid
axis([-5 5 -5 5])
hold off

% AM4 + RK4

% RK4
figure
plot(a(1,:),'co:')
hold on
plot(a(2,:),'co:')
hold on
plot(a(3,:),'co:')
hold on
plot(a(4,:),'co:')
hold on

plot(z_out)
hold on
fill(real(z_out),imag(z_out),'r')
title('Stability Region for AM4 + RK4')
xlabel('Re(z)')
ylabel('Im(z)')
grid
axis([-5 5 -5 5])
hold off
