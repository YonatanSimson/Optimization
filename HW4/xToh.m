function [h, del] = xToh(x, M, N)

del   = x(end);
a     = x(1:M+1);
T2    = aToh(N);
h     = 0.5*T2*a;
