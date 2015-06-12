function [x, Cost] = ProjectedNewton(H, b, lb, ub, x0, maxIter, tol, tolkkt)
% Solves the following problem via Projected Newton:
%
%   minimize     (1/2)*x'*H*x + b'*x
%   subject to   lb <= x <= ub

%OUTPUT:
%   x - Constrained optimal solution
%   active - Set of active constraints
%   Cost - vector of decreasing cost
%PARAMS
maxIterInner = 1000;% For CG - for inverting Hessian. Inner loop
tolInner = 1e-8;% For CG - for inverting Hessian. Inner loop
sigma = 0.3;% For Armijo rule
beta  = 0.1;% For Armijo rule
flag = 1;% For Armijo rule: 0 - unconstrained, 1-constrained case with box constraints
%Init
alpha_k = 0.99;

f = @(x)(0.5*x'*H*x+b'*x);
gradf = @(x)(H*x+b);

%iterate
xOld = inf(length(x0), 1);
x = x0;
dim = length(x0);
Hr = ones(dim); %#ok<NASGU>
d = zeros(dim, 1);
for k = 1:maxIter,
    if (norm(x-xOld)<tol*norm(xOld))
        break;
    end
    gradf_x = gradf(x);
    % Find active set
    IPlus = ActiveSet(x, gradf_x, lb, ub, tolkkt);
    nonActive = find(~IPlus);
    active    = find(IPlus);
    % Find descent direction
    if ( isempty(nonActive) )
        d = -gradf_x; %Projected gradient descent
    else
        % Newton step
        dd = -H(nonActive,nonActive)\gradf_x(nonActive); 
        d(nonActive) = dd;
        d(IPlus)     = -gradf_x(IPlus);
        %dr = ConjGrad(Hr, -P*gradf_x, -P*gradf_x, maxIterInner, tolInner); 

    end
    
    % a_k ~ arg min(f + a*d)
    alpha_k = ArmijoRule(f, x, f(x), gradf(x), d, sigma, beta, alpha_k, lb, ub, flag);
    %update
    xOld = x;
    x = Proj_B(x + alpha_k*d, lb, ub);
    if (mod(k, 1000)==0)
        disp(['Newton Iteration: ' num2str(k)])
    end
end
disp(['ProjectedNewton converged at iteration ' num2str(k)]);
Cost = f(x);