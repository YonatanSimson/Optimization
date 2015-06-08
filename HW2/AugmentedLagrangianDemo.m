clear; close all;
%Toy Example for Augmented Lagrangian
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
tol = 1e-10;
maxIter = 100;
tolkkt = 1e-3;

%% Augmented Lagrangian
mu0 = 1;
eta0 = 1;
beta = 1.01;
MaxIterAug = 1000;
x = AugmentedLagrangian( H, c, Aeq, beq, lb, ub, maxIter, tol, tolkkt, mu0, eta0, beta, MaxIterAug );

disp('Accuracy of Augmented Lagrangian')
norm(x - x_ref)

