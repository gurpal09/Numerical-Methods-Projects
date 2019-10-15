% Gurpal Singh
% 10/4/2017
% Math 567 Homework 3 Problem 4 w/ BC u'(1) = 3


iterations = 10;
for k = 1:iterations
N = 8 * 2^k;
h = 1/N;

% Domain
x = 0:h:1; x = x';

% Exact Solution
% u = @(x) x - 0.5;
% up = @(x) 1;
% upp = @(x) 0;

u = @(x) cos(2 * pi * x) + 2;
up = @(x) 2 * pi * sin(2 * pi *x);
upp = @(x) -4 * pi^2 * cos(2 * pi * x);

% BC Variables
sigma = 0;
beta = 3;


% # of Nodes
M = length(x);

% Construct A
e = ones(M,1);
A = spdiags([e -2*e e], -1:1, M, M);

% zero out bottom and top row
A(1,:) = 0;
A(M,:) = 0;

A(1,1) = -1;
A(1,2) = 1;
A(M,M) = -1;
A(M,M-1) = 1;

% Multiply by ( 1 / h^2 )
A = (1/(h^2)) * A;

% Add column and row of ones
A = [A ones(M,1)];
A = [A; ones(1,M+1)];
A(M+1,M+1) = 0;

% Construct F
% F = zeros(M,1); % Use this for u(x) = x debug case
F = upp(x);
% F = upp(x);

% Apply BC
F(1,1) = sigma/h + (1/2)*F(1,1);
F(M,1) = -beta/h + (1/2)*F(M,1);

exact = u(x);

% Add zero to F
F = [F;0];

% Solve
U = A\F;

% Error
err = U(1:M) - exact; 
E(k,1) = max(abs(err));
Nval(k,1) = N;
hval(k,1) = h;
end

f1 = figure;
% axis([-0.1 1.1 -1 1.1])
title('Exact vs. Approx (Second Order Centered)')
xlabel('x (m)')
ylabel('U')
grid on
hold on
plot(x,exact,'r')
hold on
plot(x,U(1:M,1),'b--')
legend('Exact', 'Approx')

f2 = figure;
plot(log(hval),log(E),'ro')
title('Error vs. h')
xlabel('Log(h)')
ylabel('Log(Error)')
hold on


% Slope of line Using first and last point to show order of Accuracy
for k=1:iterations-1
slope(k) = (log(E(k+1)) - log(E(k)))/(log(hval(k+1)) - log(hval(k)));
end
slope

% Results
fprintf(" N \t  \t  h \t \t \t Error\n")
for i = 1:iterations
    fprintf(" %d \t %.8f \t %16.8e \n", Nval(i), hval(i), E(i))
end






