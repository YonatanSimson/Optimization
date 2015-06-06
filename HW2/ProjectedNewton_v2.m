function [x, Cost] = ProjectedNewton_v2(H, f, gradf, lb, ub, x0, maxIter, tol, tolkkt)
% Solves the following problem via Projected Newton:
%
%   minimize     f(x)
%   subject to   lb <= x <= ub

%OUTPUT:
%   x - Constrained optimal solution
%   active - Set of active constraints
%   Cost - vector of decreasing cost
%PARAMS
maxIterInner = 1000;% For CG - for inverting Hessian. Inner loop
tolInner = 1e-10;% For CG - for inverting Hessian. Inner loop
sigma = 0.2;% For Armijo rule
beta  = 0.6;% For Armijo rule
flag = 1;% For Armijo rule: 0 - unconstrained, 1-constrained case with box constraints
% epsilon = 1e-6;%For active set on the box constraints
%Init
alpha_k = 0.9;
Cost = zeros(1, maxIter);

%iterate
xOld = inf(length(x0), 1);
x = x0;
dim = length(x0);
Hr = ones(dim); %#ok<NASGU>

for k = 1:maxIter,
    if (norm(x-xOld)<tol)
        disp(['Converged at iteration ' num2str(k)]);
        Cost = Cost(1:k-1);
        break;
    end
    gradf_x = gradf(x);
    % Find active set
    IPlus = ActiveSet(x, gradf_x, lb, ub, tolkkt);
    active    = find(IPlus);
    nonActive = find(~IPlus);
    % Find descent direction
    if ( isempty(nonActive) )
        d = -gradf_x; %Projected gradient descent
    else
        % Newton step
        P  = sparse(1:dim, [nonActive; active;], ones(dim, 1));%permuation matrix
        Hr = eye(dim);
        Hr(1:numel(nonActive), 1:numel(nonActive)) = H(nonActive,nonActive);
        % dr = -Hr\(P*gradf_x);
        dr = ConjGrad(Hr, -P*gradf_x, -P*gradf_x, maxIterInner, tolInner); 
        d = P'*dr;
    end
    
    % a_k ~ arg min(f + a*d)
    alpha_k = ArmijoRule(f, x, f(x), gradf(x), d, sigma, beta, alpha_k, lb, ub, flag);
    %update
    xOld = x;
    x = Proj_B(x + alpha_k*d, lb, ub);
end
