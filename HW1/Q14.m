%Q14 - toy problem
clear; close all;
t =  [-4 -3 -2 -1 0 1 2 3 4]';
f1 = [0 0 0 0 0.5 1 1 1 1]';
f2 = [0 0.0025 0.0180 0.1192 0.5000 0.8808 0.9820 0.9975 1]';

Dx = CreateDerivativeOperators1D(length(t));

rng(123);
x_tag = f1;% + 0.05*rand(size(t));
A = [1 1 1 0 0 0 0 0 0;
     0 0 1 1 1 0 0 0 0;
     0 0 0 0 1 1 1 0 0;
     0 0 0 0 0 0 1 1 1;];
A = sparse(A); 
y = A*x_tag;
L = Dx; 


m = size(L, 1);
n = size(L, 2);

x_0 = 0.5*ones(length(x_tag), 1);
epsilon = 0.00001;
lambda = 1;
maxIter = 300;
alpha = 0.01;
x = x_0;
W = speye(m);
for k = 1:maxIter,
    x = x - alpha*(A'*(A*x - y) + lambda*L'*W*(L*x));
    gamma = L*x;
    gamma(abs(gamma)<epsilon) = epsilon;
    W = sparse(1:m,1:m,abs(1./gamma));
end

figure; plot(t, x, t, f1); xlabel('t'); ylabel('f(t)'); 
legend('\hat{f1}', 'f1')

