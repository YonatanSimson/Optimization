%Q13

%% define LPF parameters
L      = 500;
DeltaP = 0.1;
DeltaS = 0.001;
DeltaT = 0.05;
N      = 47;
N1     = 18;

%% Equiriple design of LPF-using LP
wp = 0.26*pi;
ws = 0.34*pi;
wc = 0.30*pi;
w = (0:L)*pi/L;

%% Lob Barrier Parameters
t0      = 500;
mu      = 15;
epsilon = 1e-16;

%% First FIR design problem
[A, b, c] = FirstFIRProblem(wp, ws, wc, L, N, DeltaP, DeltaS);

tic
x = LogBarrierSolver(A, b, c, t0, mu, epsilon);
disp('Time for log Barrier')
toc

tic;
options = optimoptions('linprog', 'Algorithm', 'simplex');
x_ref = linprog(c, A, b, [], [], [], [], [], options);
disp('Time for simplex(linprog)')
toc

disp('Difference between linprog and my Log-Barrier Method:')
norm(x-x_ref)

[h, del]   = xToh(x,     N);
[h_ref, ~] = xToh(x_ref, N);

figure;
stem(0:length(h)-1, h);
hold on;
stem(0:length(h_ref)-1, h_ref, 'g');
hold off
legend('Log-Barrier Method', 'Simplex(linprog)')
title('Impulse responce - First problem')
ylabel('h[n]');
xlabel('n');

%% Second FIR design problem with time constraints
[A, b, f] = SecondFIRProblem(wp, ws, wc, L, N, N1, DeltaP, DeltaS, DeltaT);


tic
x = LogBarrierSolver(A, b, c, t0, mu, epsilon);
disp('Time for log Barrier')
toc

tic;
x_ref = linprog(c, A, b, [], [], [], [], [], options);
disp('Time for simplex(linprog)')
toc

disp('Difference between linprog and my Log-Barrier Method:')
norm(x-x_ref)

[h, del]   = xToh(x,     N);
[h_ref, ~] = xToh(x_ref, N);

figure;
stem(0:length(h)-1, h);
hold on;
stem(0:length(h_ref)-1, h_ref, 'g');
hold off
legend('Log-Barrier Method', 'Simplex(linprog)')
title('Impulse responce - With time constraints')
ylabel('h[n]');
xlabel('n');




