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
w = (0:L)*pi/L;

%% First FIR design problem
[A, b, c, M] = FirstFIRProblem(wp, ws, wc, L, N, DeltaP, DeltaS);

%% Lob Barrier Parameters
t0 = 500;
mu = 15;
epsilon = 1e-6;
options = optimoptions('linprog', 'Algorithm', 'simplex', 'Display', 'iter');

%% Phase 1 - Find strictly feasible starting point
gamma_0 = max(-b)+ 1;
AA = [A -ones(size(A,1), 1);
      zeros(1, size(A,2)) -1;];%lower bound constraint on gamma -> -gamma<=1 -> gamma > -1
  
bb = [b; 1];
cc = [zeros(size(A,2), 1); 1];
xx_0 = [zeros(size(A,2), 1); gamma_0];

xx_feas_ref = linprog(cc, AA, bb, [], [], [], [], xx_0, options);
xx_feas = LogBarrierMethod(AA, bb, cc, t0, xx_0, mu, epsilon);

%% Phase 2
x0 = xx_feas(1:end-1);
x = LogBarrierMethod(A, b, c, t0, x0, mu, epsilon);
x_ref = linprog(c, A, b, [], [], [], [], x0, options);
disp('Difference between linprog and my Log-Barrier Method:')
norm(x-x_ref)

[h, del] = xToh(x, M, N);
[h_ref, ~] = xToh(x_ref, M, N);

figure;
stem(0:length(h)-1, h);
hold on;
stem(0:length(h_ref)-1, h_ref, 'g');
hold off
legend('Log-Barrier Method', 'Simplex(linprog)')
title('Impulse responce')
ylabel('h[n]');
xlabel('n');

