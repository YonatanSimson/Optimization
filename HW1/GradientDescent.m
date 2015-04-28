function [x, cost] = GradientDescent(A, b, x_0, numOfIter, tol)
%   minimizes the cost of using x as the parameter for ||Ax-b||^2

numOfIterForLineSearch = 100;
cost = zeros(1, numOfIter);
alphaMax = 0.1;
rng(123);
x_0 = rand(size(x_0));
x = x_0;

for k = 1:numOfIter,
    [cost(k), grad] = CalcCost(x);
    if ( norm(grad) < tol )
        cost = cost(1:k);
        break;
    end
    %line search for alpha_k
    %minimize f(x_k + a_k*d_k) as a function of a_k
    alpha_k = BisectionLineSearch(@(x)CalcCost(x), 0, alphaMax, x, -grad, numOfIterForLineSearch, tol);
    x = x - alpha_k*grad;
end


%% calculate cost and gradient function

    function [J, grad] = CalcCost(x)
        e = (A*x -b);
        J = e'*e;
        grad = 2*A'*e;
    end

end


