function [f, g, H] = myfunc(x)
f = 3*x(1)^2 + 2*x(1)*x(2) + x(2)^2;    % Cost function
g(1) = 6*x(1)+2*x(2);
g(2) = 2*x(1)+2*x(2);
H = [6 2; 2 2];

