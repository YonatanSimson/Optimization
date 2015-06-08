function [ x ] = AugmentedLagrangian( H, c, Aeq, beq, lb, ub, maxIter, tol, tolkkt, mu0, eta0, beta, MaxIterAug )
% Solves the following problem via AugmentedLagrangian:
%
%   minimize     (1/2)*x'*H*x + c'*x
%   subject to   lb <= x <= ub and Aeq*x = beq

x0 = 0.5*(lb + ub);

%init
mu = mu0;
eta = eta0;
x = x0;
xOld = inf(length(x0), 1);

CostTot = zeros(1, MaxIterAug);
CostMin = zeros(1, MaxIterAug);
EqCost = zeros(1, MaxIterAug);

%iterate
for k = 1:MaxIterAug,
    Htild = H + mu * (Aeq' * Aeq);%For augmented form
    cTild = c - mu * Aeq'* beq - eta * Aeq';%For augmented form

    if ( sum((Aeq * x - beq).^2) < tol && ...
        norm(x-xOld) < tol)
        CostTot = CostTot(1:k-1);
        CostMin = CostMin(1:k-1);
        EqCost = EqCost(1:k-1);
        
        disp('Augmented Lagragian converged at iteration')
        disp(k)
        break;
    end
    [x, CostTot(k)] = ProjectedNewton(Htild, cTild, lb, ub, x, maxIter, tol, tolkkt);

%     f = @(x)(0.5* x' * H * x + c' * x + mu/2 * (Aeq * x - beq)^2  - eta * (Aeq * x - beq));
%     gradf = @(x)(H * x + c + mu * ((Aeq' * Aeq) * x - Aeq'* beq) - eta * Aeq') ;
%     [x, CostTot(k)] = GradientProjection(f, gradf, lb, ub, x, maxIter, tol);
%     [x, CostTot(k)] = quadprog(Htild, cTild, ...
%                      [], [], ...
%                      [], [], ...%equality cond
%                      lb, ub, ...%box constraints
%                      x, ...%starting point
%                      options) ;%instead of projected Newton
    EqCost(k) = sum((Aeq * x - beq).^2);
    CostMin(k) = 0.5* x' * H * x + c' * x;
    %update mu,eta
    eta = eta - mu*(Aeq * x - beq);
    mu  = mu*beta;

    disp('Augmented Lagragian iteration')
    disp(k)
    xOld = x;

end


end

