%%TestNewton

%Toy Example for Armijo rule
H=[5 -2 -1; -2 4 3; -1 3 5];
c=[2; -35; -47];
lb = [6; 6; 6];
ub = [8; 8; 8];
options = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','final-detailed');
x_unc = quadprog(H,c,[],[],[], [], [],[], [], options);

maxIterInner = 1000;
tolInner = 1e-10;
tolkkt = 1e-3;
tol = 1e-8;
maxIter = 1000;
x_cg = ConjGrad(H,-c,[3; 3; 3], maxIter, tol);

%% my Newton Step
f = @(x)(0.5*x'*H*x+c'*x);
grad_f = @(x)(H*x+c);
hessian_f = @(x)(H);
x0 = [0; 0; 0];
x = NewtonMethod_v2(f, grad_f, hessian_f, x0, maxIter, tol);