function [A, b, c, M] = FirstFIRProblem(wp, ws, wc, L, N, DeltaP, DeltaS)

w = (0:L)*pi/L;
%weighting matrix
S = diag(1/DeltaP*(w<=wp) + 1/DeltaS*(w>=ws));
d = (w<=wc)';%desired LPF frequency responce


M = (N-1)/2;
C = [ones(L+1, 1) cos(w'*(1:M))];

A = [ S*C -ones(L+1,1);
     -S*C -ones(L+1,1);];
b = [S*d; -S*d];
c = [zeros(M+1,1); 1];