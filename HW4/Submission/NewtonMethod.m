function [x, Cost] = NewtonMethod(A, b, c, t, x0, maxIter, tol)
% Solves the following centering problem via Projected Newton:
%
%   minimize     t*c^T*x + phi(x)
%   subject to   Ax < b

%OUTPUT:
%   x - Constrained optimal solution
%   Cost - vector of decreasing cost
%PARAMS
sigma        = 0.01;% For Armijo rule
beta         = 0.5;% For Armijo rule
maxIterInner = 10000;% For CG - for inverting Hessian. Inner loop
tolInner     = 1e-6;% For CG - for inverting Hessian. Inner loop

%Init
alpha_0      = 1;% has to be 1 for Newton method to work

f         = @(x)(t*c'*x -sum(log(b - A * x)));
grad_f    = @(x)(t*c + sum(bsxfun(@rdivide, A', (b - A * x)'), 2));
hessian_f = @(x)(A' * diag(1./((b - A * x).^2))*A);

%iterate
CostOld = inf;
xOld    = inf(length(x0), 1);
x       = x0;
dim     = length(x0);
H       = ones(dim); %#ok<NASGU>
d       = inf(dim, 1); %#ok<NASGU>
Costs   = zeros(1, maxIter);
m       = size(A, 1); 
alpha_k = 1;
for k = 1:maxIter,
    Cost = f(x);
    if ( norm(x-xOld)<tol*norm(xOld) )
        break;
    end
    gradf_x = grad_f(x);
    H = hessian_f(x);
    % Newton step

    
%     d = -H\gradf_x; 
    
    %Truncated Newton is faster and works better when H is singular
    d = ConjGrad(H, -gradf_x, -gradf_x, maxIterInner, tolInner); 


    % a_k ~ arg min(f + a*d)
    alpha_k = ArmijoRule(f, A, b, x, Cost, gradf_x, d, sigma, beta, alpha_0);
    
    %update
    xOld = x;
    x = x + alpha_k*d;
    
    if ( sum(A*x < b) < m )
        error('Point returned from newton step is not strictly feasible'); 
    end

    if (mod(k, 1000)==0)
        disp(['Newton Iteration: ' num2str(k)])
    end
end
disp(['Newton step converged at iteration ' num2str(k)]);

