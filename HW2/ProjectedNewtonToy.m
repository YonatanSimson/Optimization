%Toy Example for Armijo rule

H=[5 -2 -1; -2 4 3; -1 3 5];
c=[2; -35; -47];
lb = [6; 6; 6];
ub = [8; 8; 8];
options = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','final-detailed');
x_ref = quadprog(H,c,[],[],[], [], lb,ub, [], options);
x_unc = quadprog(H,c,[],[],[], [], [],[], [], options);
maxIterInner = 1000;
tolInner = 1e-10;
x_cg = ConjGrad(H,-c,[3; 3; 3], maxIterInner, tolInner);

f = @(x)(0.5*x'*H*x+c'*x);

gradf = @(x)(H*x+c);

%params
maxIter = 1000;%For Projected Newton
tol = 1e-8;%For Projected Newton
maxIterInner = 1000;% For CG - for inverting Hessian. Inner loop
tolInner = 1e-10;% For CG - for inverting Hessian. Inner loop
sigma = 0.2;% For Armijo rule
beta  = 0.6;% For Armijo rule
flag = 1;%0 - unconstrained, 1-constrained case with box constraints
epsilon = 1e-6;%for active set
%Init
x0 = 0.5*(lb + ub);
alpha_k = 0.9;
%iterate
xOld = inf(length(x0), 1);
x = x0;
dim = length(x0);
Hr = ones(dim);

for k = 1:maxIter,
    if (norm(x-xOld)<tol)
        disp(['Converged at iteration ' num2str(k)]);
        break;
    end
    gradf_x = gradf(x);
    % Find active set
    IPlus = ActiveSet(x, gradf_x, lb, ub, 1e-6);
    active    = find(IPlus);
    nonActive = find(~IPlus);
    P  = sparse(1:dim, [nonActive; active;], ones(dim, 1));%permuation matrix
    Hr = eye(dim);
    Hr(1:numel(nonActive), 1:numel(nonActive)) = H(nonActive,nonActive);
    % Find descent direction
    if ( isempty(nonActive) )
        d = -gradf_x; %Projected gradient descent
    else
        % Newton step
        dr = ConjGrad(Hr, -P*gradf_x, -P*gradf_x, maxIterInner, tolInner); % dr = -Hr\(P*gradf_x);
        d = P'*dr;
    end
    
    % a_k ~ arg min(f + a*d)
    alpha_k = ArmijoRule(f, x, f(x), gradf(x), d, sigma, beta, alpha_k, lb, ub, flag);
    %update
    xOld = x;
    x = Proj_B(x + alpha_k*d, lb, ub);
end

