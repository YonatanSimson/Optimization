%IRLS

function [x, cost] = IRLS(A, L, y, epsilon, lambda, maxIter, tol)
m = size(L, 1);

W = speye(m);
cost = zeros(3, maxIter);

maxIterIRLS = 10;
%Solution using lsqr
B = [y; zeros(size(L, 1), 1)];
for k = 1:maxIterIRLS,
    AA = [A; sqrt(lambda)*sqrt(W)*L];
    x = lsqr(AA, B, tol, maxIter);%replace with CG

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
        val = 0.5*(e'*e) + lambda*sum(abs(L*x));
    end

end


