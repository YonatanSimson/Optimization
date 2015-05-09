%IRLS

function [x, cost] = IRLS(A, L, y, epsilon, alpha, maxIter, tol)
m = size(L, 1);
n = size(L, 2);

W = speye(m);
cost = zeros(3, maxIter);

maxIterIRLS = 10;
B = [y; zeros(size(L, 1), 1)];
rng(123);
x = 0.5*rand(n, 1);

for k = 1:maxIterIRLS,
    AA = [A; sqrt(alpha)*sqrt(W)*L];
%     x_tilde = lsqr(AA, B, tol, maxIter);%replace with CG
    x = CG_LS1(x,AA,B,tol,maxIter);

    cost(1,k) = f(x);
    cost(2,k) = norm(AA*x - B);
    cost(3,k) = norm(L*x, 1);
    if (k > 1 && cost(1, k) > cost(1, k-1))
        error('IRLS diverges')
    end
    
    gamma = abs(L*x);
    gamma(gamma<epsilon) = epsilon;
    W = sparse(1:m,1:m,1./gamma);
end


%% Nested function - for cost evaluation
    function val = f(x)
        e = (AA*x - B);
        val = 0.5*(e'*e) + alpha*sum(abs(L*x));
    end

end


