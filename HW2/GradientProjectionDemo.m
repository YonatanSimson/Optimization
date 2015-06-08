clear; close all;
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
tolkkt = 1e-3;
tol = 1e-8;
maxIter = 100;
x_cg = ConjGrad(H,-c,[3; 3; 3], maxIter, tol);

%% my Projected Newton
x0 = 0.5*(lb + ub);
f = @(x)(0.5*x'*H*x+c'*x);
gradf = @(x)(H*x+c);

[x_gp, Cost] = GradientProjection(f, gradf, lb, ub, x0, maxIter, tol);
disp('Accuracy of Gradient Projection')
norm(x_gp - x_ref)

