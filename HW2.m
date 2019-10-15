E = zeros(10,1);
Nval = zeros(10,1);
scaled_c_num = zeros(10,1);
for k = 1:10
% Parameters
n = k;
L = 1;
N = 8 * 2^(n);
h = 1/N;

% Domain
x = 0:h:1; x = x';

% # of Nodes
M = length(x);

% Construct A
e = ones(M,1);
A = spdiags([e -2*e e], -1:1, M, M);
A = (1 / (h^2)) * A;

% U & F
U = zeros(M,1);
F = zeros(M,1);

% Boundary Conditions
U(1,1) = 0;
U(M,1) = 0;

% Evaluating f(x)
F = -4 * pi^2 * sin(2 * pi * x);

% Solving for Approximate Solution
U_New = A\F;

% Solving for Exact Solution
exact = zeros(M,1);
exact = sin(2 * pi * x); 

 
 Nval(k,1) = N;
 infA = norm(inv(A),inf);
 tao = - (A * (U_New - exact)) ;
 E(k,1) = norm(inv(A) * tao, inf);
 
% Part 2
% Condition Number
c_num = condest(A);
scaled_c_num(k) = eps(1) * c_num;
 
end



% figure
% title('Exact Solution')
% xlabel('x')
% ylabel('u(x)')
% hold on
% plot(x,exact)
% hold on
% plot(x,zeros(length(x)), 'r')



