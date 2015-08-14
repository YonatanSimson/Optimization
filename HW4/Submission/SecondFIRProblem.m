function [A, b, f] = SecondFIRProblem(wp, ws, wc, L, N, N1, DeltaP, DeltaS, DeltaT)

M = (N-1)/2;

%% Frequency
w = (0:L)*pi/L;

%weighting matrix
S = diag(1/DeltaP*(w<=wp) + 1/DeltaS*(w>=ws));
d = (w<=wc)';%desired LPF frequency responce

%% Get time constraints
T1 = fliplr(tril(ones(N)));
T2 = aToh(N);
T = T1(1:N1+1, :)*T2;

%% Calc LP pramaters A, b, f

C = [ones(L+1, 1) cos(w'*(1:M))];

A = [ S*C -ones(L+1,1);
     -S*C -ones(L+1,1);
      T   zeros(N1+1,1);
     -T   zeros(N1+1,1);];
v = 2*DeltaT*ones(N1+1, 1);
b = [S*d; -S*d; v; v];
f = [zeros(M+1,1); 1];
