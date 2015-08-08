%Q13



%% define LPF parameters
L = 500;
DeltaP = 0.1;
DeltaS = 0.001;
N = 47;

%% Equiriple design of LPF-using LP
wp = 0.26*pi;
ws = 0.34*pi;
wc = 0.30*pi;

%% First FIR design problem

c = [-5; -4; -6];
A =  [1 -1  1
      3  2  4
      3  2  0
      -1 0  0
      0  -1 0
      0  0 -1];
b = [20; 42; 30; 0; 0; 0;];

[x,fval,exitflag,output,lambda] = linprog(c,A,b,[],[],[]);

% [A, b, c, M] = FirstFIRProblem(wp, ws, wc, L, N, DeltaP, DeltaS);

%% Lob Barrier Parameters
t0 = 500;
mu = 10;
epsilon = 1e-6;
%% Phase 1
gamma_0 = max(b)+ eps;
AA = [A -ones(size(A,1), 1);
      zeros(1, size(A,2)) -1;];%lower bound constraint on gamma -> -gamma<=1
  
bb = [b; 1];
cc = [zeros(size(A,2), 1); 1];
xx_0 = [zeros(size(A,2), 1); gamma_0];
options = optimoptions('linprog', 'Algorithm', 'simplex');
xx_feas = linprog(cc, AA, bb, [], [], [], [], [], options);

xx_feas = LogBarrierMethod(AA, bb, cc, t0, xx_0, mu, epsilon);

%check to see if we really a have strictly feasible solution
x0 = xx_feas(1:end-1);
tt = A*x0 < b;

%% Phase 2
x0 = xx_feas(1:end-1);
x = LogBarrierMethod(A, b, c, t0, x0, mu, epsilon);



