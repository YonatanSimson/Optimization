function [x, CostFinal] = GradientProjection(f, gradf, lb, ub, x0, maxIter, tol)
% Solves the following problem via Projected Newton:
%
%   minimize     f(x)
%   subject to   lb <= x <= ub

%OUTPUT:
%   x - Constrained optimal solution
%   active - Set of active constraints
%   Cost - vector of decreasing cost
%PARAMS
sigma = 0.3;% For Armijo rule
beta  = 0.1;% For Armijo rule
flag = 1;% For Armijo rule: 0 - unconstrained, 1-constrained case with box constraints
%Init
alpha_k = 0.9;

%iterate
xOld = inf(length(x0), 1);
x = x0;
for k = 1:maxIter,
    d = -gradf(x);
    if (norm(x-xOld)<tol*norm(xOld))
        break;
    end
    %a_k = arg min(f + a*d)
    alpha_k = ArmijoRule(f, x, f(x), gradf(x), d, sigma, beta, alpha_k, lb, ub, flag);
    %update
    xOld = x;
    x = Proj_B(x + alpha_k*d, lb, ub);
end
disp(['GradientProjection Converged at iteration ' num2str(k)]);
CostFinal = f(x);