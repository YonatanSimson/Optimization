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
[A, b, c, M] = FirstFIRProblem(wp, ws, wc, L, N, DeltaP, DeltaS);

%% Lob Barrier Parameters
t0 = 500;
mu = 10;
epsilon = 1e-6;
%% Phase 1
gamma_0 = max(abs(b)) + 10*eps;
AA = [A -ones(size(A,1), 1)];
cc = [zeros(size(A,2), 1); 1];
xx_0 = [zeros(size(A,2), 1); gamma_0];
xx_feas = LogBarrierMethod(AA, b, cc, t0, xx_0, mu, epsilon);

%% Phase 2
x0 = xx_feas(1:end-1);
x = LogBarrierMethod(A, b, c, t0, x0, mu, epsilon);



