function [h, del] = xToh(x, N)
M     = (N-1)/2;
del   = x(end);
a     = x(1:M+1);
T2    = aToh(N);
h     = 0.5*T2*a;
