function [J, grad] = costFunction(x, A, L, y, lambda)
%COSTFUNCTION Compute cost and gradient for gradient descent
%   J = COSTFUNCTION(x, A, y) computes the cost of using x as the
%   parameter for ||Ax-y||^2 + lambda*||Lx||^2

Q = (A'*A)+lambda*(L'*L);
J = x'*Q*x - 2*y'*A*x + y'*y;
grad = 2*Q*x - 2*A'*y;

end