% c = [1 -2 1];
% x = [-1 0 1];
% p = 2;
% k = 2;
c = [-1/6, 2, -13/2, 28/3, -13,2, 2, -1/6];
x = [-3, -2, -1, 0, 1, 2, 3];
k = 4;
p = 6;

n = p + k;

% Computing Coefficient for Leading order term
sigma = 0;
for j = 1:length(x)
   sigma = sigma + (c(j)*x(j)).^n;
end

% Loop to go through n values
for i = 0:n
    
    if (i < k) 
        S = 0;
    elseif (i == k)
        S = factorial(i);
    elseif (k < n) && ((i < p+k))
        S = 0;
    elseif (i == p + k)
        S = sigma;
    end
    fprintf("n = %d \t  S_n = %d\n",i,S);
end

% Print Results
fprintf("\n");
fprintf("----------- Leading Order Term -----------\n");
fprintf("Leading Order term: (%s) * u^(%d) * h^(%d)\n", strtrim(rats(S/factorial(n))),n,p)