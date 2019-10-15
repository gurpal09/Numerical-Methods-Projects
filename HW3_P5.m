% Gurpal Singh
% 10/07/2017
% Math 567 HW 3 #5a

close all
clear all; clc
global kmax tol beta_minus beta_plus epsilon AMatrix
global Sol A1 B1 p q

beta_minus = 1;
beta_plus = 1;
epsilon = 0.5;

A1 = 1;
B1 = 1;
p = 2;
q = 2;

u = @(x,y) A1*x^p + B1*y^q;

% f = @(x,y) -100 * exp(-100*((x - 0.5)^2) + (y - 0.5)^2);
f = @(x,y) A1*(p-1)*p*x^(p-2) + B1*(q-1)*q*y^(q-2);

beta = @(x,y) beta_minus + (beta_plus - beta_minus) * ((1 + tanh((x-0.5)/epsilon)) / 2);



Nvals = [10 25 50 75 100 125 150];
iterations = length(Nvals);
taoVals = zeros(iterations,1);

for k = 1:iterations
N = Nvals(k);
xe = linspace(0,1,N+1)';
x = xe(2:end-1);
y = x;
h = 1/N;
[xm,ym] = meshgrid(x,y);

B = zeros(N-1, N-1);
B(1,:) = 1;
B(end,:) = 1;
B(:,1) = 1;
B(:,end) = 1;
B = (1 / h^2) * B;
F = f(xm,ym) - B;

tol = 5e-10;
kmax = 1000;

% % Conjugate Gradient
[U,E_cg] = conj_gradient(F);
tao = norm(AMatrix - Sol,inf);
taoVals(k,1) = tao;
hvals(k,1) = h;
end

figure
title('Error vs. h')
xlabel('Log(h)')
ylabel('Log(Error)')
hold on
plot(log(hvals),log(taoVals),'ro')

% Slope of line Using first and last point to show order of Accuracy
for k=1:iterations-1
slope(k) = (log(taoVals(k+1)) - log(taoVals(k)))/(log(hvals(k+1)) - log(hvals(k)));
end
slope

% Results
fprintf(" N \t  \t  h \t \t \t Error\n")
for i = 1:k
    fprintf(" %d \t %.8f \t %.8f \n", Nvals(i), hvals(i), taoVals(i))
end



% hd = surf(xm,ym,U);
% set(hd,'edgecolor','none');
% set(hd,'facecolor','interp');
% title('\beta = 1000')
% axis([0 1 0 1]);
% daspect([1 1 1]);
% 
% figure;
% clf;
% semilogy(E_jacobi,'b.','markersize',10);
% 
% semilogy(E_sd,'g.','markersize',10);
 
% semilogy(E_cg,'k.-','markersize',10);
 
% semilogy(E_pcg,'mo-','markersize',8);

function Au = matvec_pcg(U)

    Nm1 = sqrt(size(U,1));
    U = reshape(U,Nm1,Nm1);
    AU = matvec(U);
    
    Au = -AU(:);

end

function Au = matvec(U,Ubig)

% Thermal Conductivity
global epsilon beta_plus beta_minus A1 B1 p q
beta = @(x,y) beta_minus + (beta_plus - beta_minus) * ((1 + tanh((x-0.5)/epsilon)) / 2);

u = @(x,y) A1*x^p + B1*y^q;
% f = @(x,y) -100 * exp(-100*((x - 0.5)^2) + (y - 0.5)^2);
f = @(x,y) A1*(p-1)*p*x^(p-2) + B1*(q-1)*q*y^(q-2);

beta = @(x,y) beta_minus + (beta_plus - beta_minus) * ((1 + tanh((x-0.5)/epsilon)) / 2);

lap_u = @(x,y) A1*(p-1)*p*x^(p-2) + B1*(q-1)*q*y^(q-2);

dbeta_dx = @(x,y) ((beta_plus - beta_minus) * sech((x - 0.5)/epsilon)^2) / 2 * epsilon;
dbeta_dy = @(x,y) 0;

du_dx = @(x,y) A1 * p * x^(p-1);
du_dy = @(x,y) B1 * q * y^(q-1);

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

Test = zeros(N-1,N-1);
for i = 2:N

    for j = 2:N
        x_beta = xe(i);
        y_beta = ye(j); 
        
    %Beta Weights
    Beta_Up = beta(x_beta,y_beta);
    Beta_Bottom = beta(x_beta,y_beta);
    Beta_Left = beta(x_beta - h/2,y_beta);
    Beta_Right = beta(x_beta + h/2,y_beta);
    Beta_Center = beta(x_beta,y_beta);
         
        Au(i-1,j-1) = (Beta_Right*Ubig(i+1,j) + Beta_Left*Ubig(i-1,j) + Beta_Up*Ubig(i,j+1) + ...
            Beta_Bottom*Ubig(i,j-1) - 4*Beta_Center*Ubig(i,j))/h^2;
        
        % Exact Calculation
        Test(i-1,j-1) = beta(x_beta,y_beta) * lap_u(x_beta,y_beta) + ((dbeta_dx(x_beta,y_beta) * du_dx(x_beta,y_beta)) + (dbeta_dy(x_beta,y_beta) * du_dy(x_beta,y_beta)));
    end
    
    
end
   global AMatrix Sol
   AMatrix = Au;
   Sol = Test;
   
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
            break;
        end
        bk = dot(rkp1(:),rkp1(:))/dot(rk(:),rk(:));
        pk = rkp1 + bk*pk;
        rk = rkp1;
        uk = ukp1;
    end
    u = ukp1;
end


