function [J, grad] = costFunction(x, A, L, y, lambda)
%COSTFUNCTION Compute cost and gradient for gradient descent
%   J = COSTFUNCTION(x, A, y) computes the cost of using x as the
%   parameter for ||Ax-y||^2 + lambda*||Lx||^2

c = y'*y;
Q = (A'*A)+lambda*(L'*L);
J = x'*Q*x - 2*y'*A*x + c;
grad = 2*A'*(A*x - y) + 2*L'*(L*x);

end