clear; close all;
%Toy Example for Armijo rule
H=[5 -2 -1; -2 4 3; -1 3 5];
c=[2; -35; -47];
Aeq = [1 1 1];
beq = 19;
lb = [6; 6; 6];
ub = [8; 8; 8];
options = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','final-detailed');
x_ref = quadprog(H,c,[],[],Aeq, beq, lb,ub, [], options);
x_unc = quadprog(H,c,[],[],[], [], [],[], [], options);
tol = 1e-8;
maxIter = 100;
x_cg = ConjGrad(H,-c,[3; 3; 3], maxIter, tol);

%% Augmented Lagrangian
mu0 = 1;
eta0 = 1;
beta = 1.01;
MaxIterAug = 1000;
x0 = 0.5*(lb + ub);

%init
mu = mu0;
eta = eta0;
x = x0;
xOld = inf(length(x0), 1);

CostTot = zeros(1, MaxIterAug);
CostMin = zeros(1, MaxIterAug);
EqCost = zeros(1, MaxIterAug);
options = optimset('Algorithm', 'trust-region-reflective',...%either equality or box
    'Display','final-detailed', ...
    'LargeScale', 'on', ...
    'TolFun',tol, ...
    'MaxIter', maxIter);

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
    f = @(x)(0.5* x' * H * x + c' * x + mu/2 * (Aeq * x - beq)^2  - eta * (Aeq * x - beq));
    gradf = @(x)(H * x + c + mu * ((Aeq' * Aeq) * x - Aeq'* beq) - eta * Aeq') ;
%     [lambda, CostTot(k)] = ProjectedNewton(H, b, lb, ub, lambda, 2000, tol, tolkkt);
%     [lambda, CostTot(k)] = ProjectedNewton_v2(H, f, gradf, b, lb, ub, lambda, 2000, tol, tolkkt);

%     [x, CostTot(k)] = GradientProjection(f, gradf, lb, ub, x, maxIter, tol);
    [x, CostTot(k)] = quadprog(Htild, cTild, ...
                     [], [], ...
                     [], [], ...%equality cond
                     lb, ub, ...%box constraints
                     x, ...%starting point
                     options) ;%instead of projected Newton
    EqCost(k) = sum((Aeq * x - beq).^2);
    CostMin(k) = 0.5* x' * H * x + c' * x;
    %update mu,eta
    eta = eta - mu*(Aeq * x - beq);
    mu  = mu*beta;

    disp('Augmented Lagragian iteration')
    disp(k)
    disp('Distance from ref:')
    norm(x - x_ref)
    xOld = x;

end

disp('Accuracy of Augmented Lagrangian')
norm(x - x_ref)

