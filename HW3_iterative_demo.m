
close all

global kmax tol

u = @(x,y) cos(2*pi*x).*sin(2*pi*y);
uxx = @(x,y) -(2*pi)^2*cos(2*pi*x).*sin(2*pi*y);
uyy = @(x,y) -(2*pi)^2*cos(2*pi*x).*sin(2*pi*y);
f = @(x,y) uxx(x,y) + uyy(x,y);

N = 8;
xe = linspace(0,1,N+1)';
x = xe(2:end-1);
y = x;
h = 1/N;

[xm,ym] = meshgrid(x,y);

B = compute_bdry(u,N);

F = f(xm,ym) - B;

tol = 1e-10;
kmax = 5000;

% Jacobi
[U,E_jacobi] = splitting(F,'jacobi');


% Gauss-seidel
[U,E_gs] = splitting(F,'gs');

% Steepest-descent
[U,E_sd] = steepest_descent(F);

% Conjugate Gradient
[U,E_cg] = conj_gradient(F);

hd = surf(xm,ym,U);
set(hd,'edgecolor','none');
set(hd,'facecolor','interp');

axis([0 1 0 1]);
daspect([1 1 1]);

[U,flag,relres,iter,E_pcg] = pcg(@matvec_pcg,-F(:),tol,kmax);
fprintf('PCG flag : %d\n',flag);
fprintf('PCG iter : %d\n',iter);

Ue = u(xm,ym);
fprintf('Norm(U-Ue) : %12.4e\n',norm(U(:)-Ue(:),Inf));


figure;
clf;
semilogy(E_jacobi,'b.','markersize',10);

hold on;
semilogy(E_gs,'r.','markersize',10);
 
semilogy(E_sd,'g.','markersize',10);
 
semilogy(E_cg,'k.-','markersize',10);
 
semilogy(E_pcg,'mo-','markersize',8);

lh = legend('Jacobi','Gauss-Seidel','Steepest-descent','Conj. Grad.','PCG');
set(lh,'location','southwest')
title('Iterative Methods')


set(gca,'fontsize',16)
xlim([0 N])




function Au = matvec_pcg(U)

    Nm1 = sqrt(size(U,1));
    U = reshape(U,Nm1,Nm1);
    AU = matvec(U);
    
    Au = -AU(:);

end

function Au = matvec(U,Ubig)

N = size(U,1) + 1;
h = 1/N;

if (nargin == 1)
    Ubig = zeros(N+1,N+1);
end

Ubig(2:end-1,2:end-1) = U;

Au = zeros(N-1,N-1);
for i = 2:N
    for j = 2:N
        Au(i-1,j-1) = (Ubig(i+1,j) + Ubig(i-1,j) + Ubig(i,j+1) + ...
            Ubig(i,j-1) - 4*Ubig(i,j))/h^2;
    end
end
   
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
        fprintf('%5d %12.4e\n',k,e(k));
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