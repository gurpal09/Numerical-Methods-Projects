% Gurpal Singh
% 10/07/2017
% Math 567 HW5 5b

close all
clear all; clc
global kmax tol beta_minus beta_plus epsilon AMatrix Sol A1 B1 p q
global Total_Iterations

% Coefficients for RHS
A1 = 1;
B1 = 1;1
q = 2;


u = @(x,y) A1*x^p + B1*y^q;
f = @(x,y) -100 * exp(-100*((x - 0.5)^2) + (y - 0.5)^2);
% f = @(x,y) A1*(p-1)*p*x^(p-2) + B1*(q-1)*q*y^(q-2);

% Convergence Variables
Nval = [20 50 75 100 160];
iterations = 2;
Itval = zeros(iterations,1);

% Beta and Epsilon Dependence
Epsilon_Vals = [0.5 0.1 0.05 0.01 0.005 0];
Beta_Vals = [1 10 100 1000];

for k = 1:length(Epsilon_Vals)
    for j = 1:length(Beta_Vals)
    
        % Parameters for Beta
        beta_minus = 1;
        beta_plus = Beta_Vals(j);
        epsilon = Epsilon_Vals(k);

        N = 160;
        xe = linspace(0,1,N+1)';
        x = xe(2:end-1);
        y = x;
        h = 1/N;

        [xm,ym] = meshgrid(x,y);

        % Enforce Dirichlet BC
        B = zeros(N-1, N-1);
        B(1,:) = 1;
        B(end,:) = 1;
        B(:,1) = 1;
        B(:,end) = 1;
        B = (1 / h^2) * B;
        F = f(xm,ym) - B;

        % Set Tolerance and Max Iteration
        tol = 1e-5;
        kmax = 5000;

        % Conjugate Gradient
        [U,E_cg] = conj_gradient(F);

        Itval(k,1) = Total_Iterations;
        fprintf("Epsilon = %.4f\t Beta = %.4f\t Iterations = %d\n",epsilon, beta_plus,Total_Iterations);
        % Calculate Error
        % tao(k,1) = norm(AMatrix - Sol,inf);
      end
end

% Surface Plot
% hd = surf(xm,ym,U);
% set(hd,'edgecolor','none');
% set(hd,'facecolor','interp');
% axis([0 1 0 1]);
% daspect([1 1 1]);

% % Figure for Convergence
% figure;
% clf; 
% semilogy(E_cg,'k.-','markersize',10);

function Au = matvec_pcg(U)

    Nm1 = sqrt(size(U,1));
    U = reshape(U,Nm1,Nm1);
    AU = matvec(U);
    
    Au = -AU(:);

end

function Au = matvec(U,Ubig)

% Thermal Conductivity
global epsilon beta_plus beta_minus A1 B1 p q

current_beta = 0; % Used for piecewise epsilon = 0 case

if (epsilon > 0)
    beta = @(x,y) beta_minus + (beta_plus - beta_minus) * ((1 + tanh((x-0.5)/epsilon)) / 2);
else
    beta = @(x,y) current_beta;
end
dbeta_dy = @(x,y) 0;
dbeta_dx = @(x,y) (beta_plus - beta_minus)*sec((x-0.5)/epsilon) / (2 * epsilon);

du_dx = @(x,y) A1 * p * x^(p-1);
du_dy = @(x,y) B1 * q * y^(q-1);

laplacian_u = @(x,y) A1 * (p-1) * p * x^(p-2) + B1 * (q-1) * q * y^(q-2);

exact = @(x,y) beta(x,y) * laplacian_u(x,y) + (dbeta_dx(x,y) * du_dx(x,y) + dbeta_dy(x,y) * du_dy(x,y));

N = size(U,1) + 1;
h = 1/N;

xe = linspace(0,1,N+1)';
ye = xe;
x = xe(2:end-1);

if (nargin == 1)
    Ubig = zeros(N+1,N+1);
end

Ubig(2:end-1,2:end-1) = U;

Au = zeros(N-1,N-1);
True = zeros(N-1,N-1);
for i = 2:N
    
    
    
    for j = 2:N
        x_beta = xe(i);
        y_beta = ye(j); 
        if (epsilon == 0)
%             disp('HIII')
        if((x_beta) >= 0.5)
            current_beta = beta_plus;
        else
            current_beta = beta_minus;
        end
        end
    %Beta Weights
    Beta_Up = beta(x_beta,y_beta);
    Beta_Bottom = beta(x_beta,y_beta);
    Beta_Left = beta(x_beta - h/2,y_beta);
    Beta_Right = beta(x_beta + h/2,y_beta);
    Beta_Center = beta(x_beta,y_beta);
         
        Au(i-1,j-1) = (Beta_Right*Ubig(i+1,j) + Beta_Left*Ubig(i-1,j) + Beta_Up*Ubig(i,j+1) + ...
            Beta_Bottom*Ubig(i,j-1) - 4*Beta_Center*Ubig(i,j))/h^2;
        
       % True Solution
       True(i-1,j-1) = exact(x_beta,y_beta);
    end
    
    
end
   global AMatrix Sol
   AMatrix = Au;
   Sol = True;
end


function B = compute_bdry(u,N)

xe = linspace(0,1,N+1);
ye = xe;
[xem,yem] = meshgrid(xe,ye);

Ue = u(xem,yem);

U = zeros(N-1,N-1);

B = matvec(U,Ue);

end

function Z = scale_by_M(R,m)

N = size(R,1)+1;
h = 1/N;

switch m
    case 'none'
        Z = R;
    case 'jacobi'
        d = -4/h^2;
        Z = R/d;
    case 'gs'
        Zbig = zeros(N+1,N+1);
        Zbig(2:end-1,2:end-1) = R;
        for i = 2:N
            for j = 2:N
                Zbig(i,j) = (Zbig(i-1,j) + Zbig(i,j-1) - h^2*R(i-1,j-1))/4;
            end
        end
        Z = Zbig(2:end-1,2:end-1);
end

end


function [u,e] = splitting(F,m)

global kmax tol

uk = zeros(size(F));
rk = F - matvec(uk);
zk = scale_by_M(rk,m);

e = zeros(kmax,1);
for k = 1:kmax
    ukp1 = uk + zk;
    
    e(k) = norm(zk(:),Inf);
%     fprintf('%5d %12.4e\n',k,e(k));
    if (e(k) < tol)
        break;
    end
    rk = F - matvec(ukp1);   % Matrix vector multiply!
    zk = scale_by_M(rk,m);
    uk = ukp1;   
end
u = ukp1;

end


function [u,e] = steepest_descent(F)

    global kmax tol

    uk = zeros(size(F));
    rk = F - matvec(uk);
    
    for k = 1:kmax
        wk = matvec(rk);
        alpha = dot(rk(:),rk(:))/dot(rk(:),wk(:));
        ukp1 = uk + alpha*rk;
        rk = rk - alpha*wk;
        e(k) = norm(rk(:));
        if (e(k) < tol)
            break;
        end
        uk = ukp1;
    end
    u = ukp1;
end

function [u,e] = conj_gradient(F)

    global kmax tol
    global Total_Iterations
    uk = zeros(size(F));
    rk = F - matvec(uk);
    pk = rk;
    
    for k = 1:kmax
        wk = matvec(pk);
        alpha = dot(rk(:),rk(:))/dot(pk(:),wk(:));
        ukp1 = uk + alpha*pk;
        rkp1 = rk - alpha*wk;
        e(k) = norm(rk(:));
%         fprintf('%5d %12.4e\n',k,e(k));
        
        if (e(k) < tol)
            Total_Iterations = k;
            break;
        end
        bk = dot(rkp1(:),rkp1(:))/dot(rk(:),rk(:));
        pk = rkp1 + bk*pk;
        rk = rkp1;
        uk = ukp1;
    end
    u = ukp1;
    
end


