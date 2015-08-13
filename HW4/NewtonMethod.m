function [x, Cost] = NewtonMethod(A, b, c, t, x0, maxIter, tol)
% Solves the following centering problem via Projected Newton:
%
%   minimize     t*c^T*x + phi(x)
%   subject to   Ax < b

%OUTPUT:
%   x - Constrained optimal solution
%   Cost - vector of decreasing cost
%PARAMS
maxIterInner = 1000;% For CG - for inverting Hessian. Inner loop
tolInner = 1e-4;% For CG - for inverting Hessian. Inner loop
sigma = 0.01;% For Armijo rule
beta  = 0.5;% For Armijo rule
%Init
alpha_0 = 1;

f = @(x)(t*c'*x -sum(log(b - A * x)));
grad_f = @(x)(t*c + sum(bsxfun(@rdivide, A', (b - A * x)'), 2));
hessian_f = @(x)(A' * diag(1./((b - A * x).^2))*A);

%iterate
xOld = inf(length(x0), 1);
x = x0;
dim = length(x0);
H = ones(dim); %#ok<NASGU>
d = zeros(dim, 1); %#ok<NASGU>
Costs = zeros(1, maxIter);
for k = 1:maxIter,
    if ( norm(x-xOld)<tol*norm(xOld) )
        Costs = Costs(1:k-1);
        break;
    end
    gradf_x = grad_f(x);
    H = hessian_f(x);
    % Newton step
    d = -H\gradf_x; 
    %d = ConjGrad(H, -gradf_x, -gradf_x, maxIterInner, tolInner); 

    % a_k ~ arg min(f + a*d)
    alpha_k = ArmijoRule(f, A, b, x, f(x), grad_f(x), d, sigma, beta, alpha_0);
    %update
    xOld = x;
    x = x + alpha_k*d;
    
    m     = size(A, 1); 
    if ( sum(A*x < b) < m )
        error('Point returned from newton step is not strictly feasible'); 
    end

    if (mod(k, 1000)==0)
        disp(['Newton Iteration: ' num2str(k)])
    end
    Costs(k) = f(x);
end
disp(['Newton step converged at iteration ' num2str(k)]);
Cost = f(x);
figure(1);
plot(Costs);
title('Cost');
