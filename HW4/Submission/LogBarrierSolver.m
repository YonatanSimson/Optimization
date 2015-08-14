function [x, Cost] = LogBarrierSolver(A, b, c, t0, mu, epsilon)


%% Phase 1 - Find strictly feasible starting point

gamma_0 = max(-b)+ 1;
AA = [A -ones(size(A,1), 1);
      zeros(1, size(A,2)) -1;];%lower bound constraint on gamma -> -gamma<=1 -> gamma > -1
  
bb = [b; 1];
cc = [zeros(size(A,2), 1); 1];
xx_0 = [zeros(size(A,2), 1); gamma_0];

xx_feas = LogBarrierMethod(AA, bb, cc, t0, xx_0, mu, epsilon);

% For debug
% options = optimoptions('linprog', 'Algorithm', 'simplex');
% xx_feas_ref = linprog(cc, AA, bb, [], [], [], [], xx_0, options);

%% Phase 2
x0 = xx_feas(1:end-1);
[x, Cost] = LogBarrierMethod(A, b, c, t0, x0, mu, epsilon);
